"""
snowline_physics.py — Física de snowlines: v_frag(T) y posiciones rsnow
========================================================================

Mixin para WaterworldPipeline que agrupa:
  - setup_physics()         → updater de v_frag(T) multi-especie
  - add_snowline_fields()   → fields rsnow_{sp} en el grid (guardados en HDF5)

Referencias de velocidades de fragmentación:
  - Silicatos (baseline)  : Birnstiel et al. (2012) A&A 539, A148  [v=100 cm/s]
  - H2O ice               : Gundlach & Blum (2015) ApJ 798, 34     [v=1000 cm/s]
  - CO2 ice               : Gundlach et al. (2018) MNRAS 479, 1273 [v=500 cm/s]
  - CO  ice               : Dominik & Tielens (1997) ApJ 480, 647  [v=300 cm/s]
"""

import numpy as np
import dustpy.constants as c


class SnowlinePhysicsMixin:
    """
    Mixin de física de snowlines.
    Asume que self.sim ya tiene los componentes cargados
    (add_volatile_components llamado).
    """

    # ══════════════════════════════════════════════════════════════════════════
    # Velocidad de fragmentación dependiente de temperatura
    # ══════════════════════════════════════════════════════════════════════════

    def setup_physics(self):
        """
        Construye el updater de v_frag(T) dinámicamente a partir de
        self.active_species y self.vfrag_params.

        Perfil de escalones (de mayor a menor T_sub):
          vf = vfrag_silicates                   → baseline refractario
          vf = where(T < T_H2O, v_H2O, vf)      → sobreescribe zona fría
          vf = where(T < T_CO2, v_CO2, vf)       → sobreescribe zona más fría
          vf = where(T < T_CO,  v_CO,  vf)       → sobreescribe la más fría

        Solo las especies en active_species Y en vfrag_params contribuyen
        un escalón. El resto del disco usa el baseline de silicatos.

        Ejemplos:
          active_species=["H2O"]            → un escalón   en T < 150 K
          active_species=["H2O", "CO2"]     → dos escalones en T < 150 / 70 K
          active_species=["H2O", "CO2","CO"]→ perfil completo de tres escalones

        DEBE llamarse después de add_volatile_components().
        """
        print("Configurando v_frag(T) dinámica (updater de hielo)...")

        # Filtrar: solo especies que tengan entrada en vfrag_params
        ice_species = [sp for sp in self.active_species if sp in self.vfrag_params]

        # Ordenar por T_sub DESCENDENTE para aplicar correctly los escalones
        ice_species_sorted = sorted(
            ice_species,
            key=lambda sp: self.vfrag_params[sp][0],
            reverse=True,
        )

        if ice_species_sorted:
            for sp in ice_species_sorted:
                Tsub, vf = self.vfrag_params[sp]
                print(f"  → {sp}: T < {Tsub:.0f} K  →  v_frag = {vf:.0f} cm/s  "
                      f"({vf / 100:.0f} m/s)")
        else:
            print(f"  → Sin hielo activo: v_frag = {self.vfrag_silicates:.0f} cm/s "
                  f"(silicatos) en todo el disco")

        # Capturar parámetros en el closure (evaluados una única vez aquí)
        _ice  = [(self.vfrag_params[sp][0], self.vfrag_params[sp][1])
                 for sp in ice_species_sorted]
        _base = self.vfrag_silicates

        def v_frag_variable(sim):
            """
            Perfil v_frag(T) construido dinámicamente para cada timestep.
            Baseline: silicatos (_base cm/s).
            Cada especie helada añade un escalón: donde T < T_sub → v_ice.
            """
            T  = sim.gas.T
            vf = np.full_like(T, _base)
            for T_sub, v_ice in _ice:
                vf = np.where(T < T_sub, v_ice, vf)
            return vf

        self.sim.dust.v.frag.updater.updater = v_frag_variable
        self.sim.dust.v.frag.update()

    # ══════════════════════════════════════════════════════════════════════════
    # Fields de posición de snowline (rsnow)
    # ══════════════════════════════════════════════════════════════════════════

    def add_snowline_fields(self):
        """
        Agrega un field rsnow_{sp} por cada especie volátil activa.

        Cada field devuelve el primer radio donde T < T_sub (borde interior
        de la zona fría). El valor se guarda en cada snapshot HDF5 y tras
        read.all() tiene shape (Nt,).

        Nota sobre T_bind vs T_sub:
          comp.gas.pars.Tsub almacena la ENERGÍA de enlace [K] (e.g. 5800 K
          para H2O), NO la temperatura de sublimación física (150 K).
          Por eso este método usa sus propias T_sub_K hardcodeadas, que son
          los valores observacionales estándar.

        DEBE llamarse después de add_volatile_components() y setup_physics().
        """
        print("Registrando fields de posición de snowline (rsnow_{sp})...")

        # Temperaturas de sublimación estándar por especie
        _T_sub_known = {
            "H2O": 150.0,
            "CO2":  70.0,
            "CO":   25.0,
        }

        # Solo para especies activas con T_sub conocida (evita snowlines fantasma)
        snowline_species = {
            sp: T for sp, T in _T_sub_known.items()
            if sp in self.active_species
        }

        if not snowline_species:
            print("  → Ninguna especie activa tiene snowline definida. Omitiendo.")
            return

        print(f"  → Snowlines a registrar: {list(snowline_species.keys())} "
              f"(active_species={self.active_species})")

        def make_rsnow_updater(T_sub_K):
            """
            Closure que captura T_sub_K para un updater independiente por especie.

            Algoritmo:
              - Si T < T_sub en todo el grid (disco muy frío)  → ri[0]
              - Si T >= T_sub en todo el grid (disco muy cálido) → ri[-1]
              - Si no → primer índice donde T < T_sub
            """
            def _updater(sim):
                cold_mask = sim.gas.T < T_sub_K
                if not cold_mask.any():
                    return float(sim.grid.ri[-1])   # snowline más allá del borde
                return float(sim.grid.ri[int(np.argmax(cold_mask))])
            return _updater

        for spec_name, T_sub in snowline_species.items():
            if spec_name == "H2O":
                field_name = "r_ice"
            else:
                field_name = f"rsnow_{spec_name}"
            
            self.sim.dust.addfield(
                field_name,
                0.0,
                updater=make_rsnow_updater(T_sub),
                description=f"Posición del snowline de {spec_name} [cm]",
                save=True,
            )
            getattr(self.sim.dust, field_name).update()
            r_au = float(getattr(self.sim.dust, field_name)) / c.au
            print(f"  → {field_name} (T_sub = {T_sub:.0f} K):  {r_au:.2f} AU")

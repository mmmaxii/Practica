"""
disk_chemistry.py — Química del disco: volátiles, silicatos y fields de Σ
===========================================================================

Mixin para WaterworldPipeline que agrupa:
  - add_volatile_components()   → inyecta trazadores de volátiles desde chem.txt
                                   y el componente de fondo de silicatos
  - add_ice_sigma_fields()      → registra y actualiza SigmaDust/SigmaGas por
                                   componente en el HDF5 (via diastole)

Referencias:
  - Gundlach & Blum (2015) ApJ 798, 34          — v_frag H2O ice
  - Gundlach et al. (2018) MNRAS 479, 1273       — v_frag CO2 ice
  - Birnstiel et al. (2012) A&A 539, A148        — v_frag silicatos
"""

import numpy as np
import dustpy.constants as c


class DiskChemistryMixin:
    """
    Mixin de química del disco.
    Asume que self.sim ya está inicializado (initialize_simulation llamado).
    """

    # ══════════════════════════════════════════════════════════════════════════
    # Componentes volátiles y silicatos
    # ══════════════════════════════════════════════════════════════════════════

    def add_volatile_components(self):
        """
        Inyecta los trazadores térmicos de volátiles a partir de chem.txt
        y el componente de fondo de silicatos refractarios.

        Patrón de tripodpy:
          - Componentes solo-gas  (nu_des <= 0): gas_active=True,  dust_active=False
          - Componentes híbridos  (nu_des >  0): gas_active=True,  dust_active=True
              → fijar gas.pars.nu, Tsub, mu y dust.pars.rhos DESPUÉS de addcomponent
          - Silicatos de fondo: dust_active=True, gas_active=False,  rhos=3.5 g/cm³
          - Llamar dust.rhos.update() al final para que tripodpy recalcule rhos
        """
        print("Agregando componentes (volátiles + silicatos) desde chem.txt...")

        data = np.genfromtxt(
            'chem.txt',
            dtype=None,
            names=True,
            encoding='utf-8',
            usecols=(0, 1, 2, 3, 4, 5),
            comments='#',
        )

        # ── Gas residual (H₂ + He): se va descontando especie a especie ───────
        Sig_residual = self.sim.gas.Sigma.copy()
        frac_h = 0.9118   # fracción de masa en hidrógeno
        self.sim.gas.mu = np.ones_like(self.sim.gas.mu) * (
            frac_h * c.m_p + (1 - frac_h) * 4 * c.m_p
        ) / (0.5 * frac_h + (1 - frac_h))

        print(f"  → Especies activas: {self.active_species}")

        for element in data:
            spec_name = element['Species']
            if spec_name not in self.active_species:
                continue

            # Fracción de masa relativa a la masa atómica media del disco
            mass_frac = element['Abundance'] * element['mu'] / (frac_h + 4 * (1 - frac_h))

            if element["nu_des"] <= 0:
                # ── Trazador solo-gas ─────────────────────────────────────────
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma * mass_frac,
                    element['mu'] * c.m_p,
                    dust_active=False,
                    gas_active=True,
                )

            else:
                # ── Componente híbrido (condensación/evaporación) ─────────────
                # Inicializar polvo solo donde Σ_dust > piso numérico
                factor = np.where(
                    self.sim.dust.Sigma[:, 0] > self.sim.dust.SigmaFloor[:, 0],
                    1.0, 0.0,
                )
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma * mass_frac * factor,
                    element['mu'] * c.m_p,
                    dust_active=True,
                    gas_active=True,
                    rhos=1.,
                    dust_value=self.sim.dust.SigmaFloor.copy(),
                )
                # Parámetros de condensación/evaporación — DESPUÉS de addcomponent
                comp = getattr(self.sim.components, spec_name)
                comp.gas.pars.nu    = element['nu_des']       # frecuencia de desorción [Hz]
                comp.gas.pars.Tsub  = element['T_bind']       # energía de enlace [K]
                comp.gas.pars.mu    = element['mu'] * c.m_p   # masa molecular [g]
                comp.dust.pars.rhos = 1.0                     # densidad bulk polvo [g/cm³]

            # Descuento único del gas residual (evita doble sustracción)
            Sig_residual -= self.sim.gas.Sigma * mass_frac

        # ── Silicatos refractarios (componente de fondo) ──────────────────────
        self.sim.addcomponent(
            "silicates",
            self.sim.gas.SigmaFloor,
            1.,
            dust_active=True,
            gas_active=False,
            rhos=3.5,
            dust_value=self.sim.dust.Sigma.copy(),
        )

        # ── Gas residual → componente Default (H₂) ────────────────────────────
        self.sim.components.Default.gas.Sigma   = Sig_residual * 2   # TODO: revisar factor 2
        self.sim.components.Default.gas.pars.mu = 2 * c.m_p

        # ── Recalcular densidad interna del polvo tras añadir todos ──────────
        # Necesario porque addcomponent no actualiza dust.rhos automáticamente.
        self.sim.dust.rhos.update()
        self.sim.update()

    # ══════════════════════════════════════════════════════════════════════════
    # Fields de densidad superficial por componente (HDF5)
    # ══════════════════════════════════════════════════════════════════════════

    def add_ice_sigma_fields(self):
        """
        Registra y actualiza fields de Σ_dust / Σ_gas por componente en el HDF5.

        Problema con updaters en addfield:
          Los fields del grid con updaters se resuelven antes de que
          sim.components tenga los valores del timestep actual → siempre
          capturan el estado INICIAL.

        Solución — diastole en sim.dust.updater.diastole:
          El diastole se dispara DESPUÉS de que dust (y components) termina de
          calcular el timestep. Usamos asignación in-place field[:] = value.

        Fields creados:
          grid/SigmaDust_{sp}   [g/cm²]  — para volátiles Y silicatos
          grid/SigmaGas_{sp}    [g/cm²]  — solo para volátiles (silicatos son
                                           refractarios, no tienen fase gas)

        Después de read.all() → shape (Nt, Nr) automáticamente.

        DEBE llamarse después de add_volatile_components().
        """
        print("Registrando fields SigmaDust / SigmaGas por componente...")

        Nr       = len(self.sim.grid.r)
        _all_comp = list(self.active_species) + ["silicates"]

        # ── Registrar fields (sin updater; se escriben vía diastole) ─────────
        for sp_name in _all_comp:
            self.sim.grid.addfield(
                f"SigmaDust_{sp_name}",
                np.zeros(Nr),
                updater=None,
                description=f"Σ_dust — componente {sp_name} [g/cm²]",
                save=True,
            )
        # SigmaGas solo para volátiles
        for sp_name in self.active_species:
            self.sim.grid.addfield(
                f"SigmaGas_{sp_name}",
                np.zeros(Nr),
                updater=None,
                description=f"Σ_gas — componente {sp_name} [g/cm²]",
                save=True,
            )

        # ── Diastole: escribe los valores tras cada update de dust ────────────
        _active    = list(self.active_species)   # capturado en el closure

        def _update_sigma_comp(sim):
            """
            Escribe Σ_dust y Σ_gas de cada componente en los fields del grid.
            Se dispara automáticamente después de sim.dust.update().

            Aplica SigmaFloor para evitar que valores sub-numéricos del solver
            contaminen los cálculos de composición en PebbleAccretion3.
            """
            dust_floor = sim.dust.SigmaFloor.sum(-1)   # (Nr,)
            gas_floor  = sim.gas.SigmaFloor             # (Nr,)

            # Volátiles — polvo + gas
            for sp_name in _active:
                comp = getattr(sim.components, sp_name)
                # Polvo: (Nr, Nbins) → suma bins → clamp al piso
                getattr(sim.grid, f"SigmaDust_{sp_name}")[:] = np.maximum(
                    comp.dust.Sigma.sum(-1), dust_floor
                )
                # Gas: (Nr,) → clamp al piso
                getattr(sim.grid, f"SigmaGas_{sp_name}")[:] = np.maximum(
                    comp.gas.Sigma, gas_floor
                )

            # Silicatos — solo polvo (refractarios, sin fase gas)
            sil_val = sim.components.silicates.dust.Sigma.sum(-1)
            sim.grid.SigmaDust_silicates[:] = np.maximum(sil_val, dust_floor)

        self.sim.dust.updater.diastole = _update_sigma_comp

        # Inicializar con el estado actual (t=0)
        _update_sigma_comp(self.sim)

        for sp_name in _all_comp:
            sig_d = float(getattr(self.sim.grid, f"SigmaDust_{sp_name}").max())
            print(f"  → SigmaDust_{sp_name}: max = {sig_d:.3e} g/cm²")
        for sp_name in self.active_species:
            sig_g = float(getattr(self.sim.grid, f"SigmaGas_{sp_name}").max())
            print(f"  → SigmaGas_{sp_name}:  max = {sig_g:.3e} g/cm²")

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
from .oka_interpolation import r_snow_time_cgs


class SnowlinePhysicsMixin:
    """
    Mixin de física de snowlines.
    Implementa el snowline global dinámico (Oka et al. 2011) y la
    asignación de v_frag según la posición radial.
    """

    # ══════════════════════════════════════════════════════════════════════════
    # Velocidad de fragmentación (v_frag) basada en r_snow
    # ══════════════════════════════════════════════════════════════════════════

    def setup_physics(self):
        """
        Construye el updater de v_frag dinámicamente.
        - Interior al snowline (r < r_snow) → material refractario seco (100 cm/s)
        - Exterior al snowline (r >= r_snow) → hielo (1000 cm/s)
        """
        print("Configurando v_frag dinámica basada en r_snow(t)...")

        def v_frag_variable(sim):
            """
            Perfil v_frag(r) actualizado dinámicamente para cada timestep.
            Depende exclusivamente de la posición global de la snowline.

            Physical interpretation of the snowline evolution:

            At very early times, the accretion rate exceeds ~1e-7 Msun/yr,
            placing the disk in an extreme viscous-heating regime. In this phase,
            the computed snowline may lie at unrealistically large radii and is
            not interpreted as physically meaningful for planet formation.

            Therefore, the snowline evolution is only considered once the disk
            accretion rate decreases below the critical threshold:

                M_dot <= 1e-7 Msun/yr

            At this stage, the snowline enters the dynamically relevant region
            of the disk with an initial radius approximately equal to the maximum
            physical ice-line radius:

                R_ice,max ~ 2.73 AU

            From that point onward, the snowline migrates inward secularly
            following the thermal evolution of the disk.
            """
            # Forzar actualización de r_snow para el timestep actual
            sim.dust.r_snow.update()
            
            # Leer el valor precalculado para este step
            r_snow = float(sim.dust.r_snow)
            
            # v_frag_silicates = 100 cm/s (usualmente)
            # v_ice ahora dependerá del parámetro v_frag_m_s que se pasa al pipeline
            v_refractory = self.vfrag_silicates
            v_ice = self.v_frag_m_s * 100.0
            
            # Condición: si r < r_snow -> polvo seco frágil, si r >= r_snow -> hielo resistente
            vf = np.where(sim.grid.r < r_snow, v_refractory, v_ice)
            return vf

        self.sim.dust.v.frag.updater.updater = v_frag_variable
        # IMPORTANTE: al usar update(), se forzará el recálculo
        self.sim.dust.v.frag.update()

    # ══════════════════════════════════════════════════════════════════════════
    # Field global de posición de snowline (r_snow)
    # ══════════════════════════════════════════════════════════════════════════

    def add_snowline_fields(self):
        """
        Agrega la propiedad global dinámica `r_snow` a la componente de polvo.
        Esta propiedad interpola los datos de Oka et al. 2011 según sim.t
        y se guarda en los HDF5.
        """
        print("Registrando field de posición de snowline dinámico (r_snow)...")

        def make_rsnow_updater():
            def _updater(sim):
                # Calcular r_snow en CGS basándonos en sim.t
                return float(r_snow_time_cgs(float(sim.t)))
            return _updater

        self.sim.dust.addfield(
            "r_snow",
            0.0,
            updater=make_rsnow_updater(),
            description="Posición dinámica del snowline de H2O [cm] (Oka et al. 2011)",
            save=True,
        )
        # Actualizar su valor inicial en t=0
        self.sim.dust.r_snow.update()
        
        r_au = float(self.sim.dust.r_snow) / c.au
        print(f"  → Snowline inicial (t=0): {r_au:.2f} AU")

"""
pipeline_snowlines.py — Pipeline principal del modelo PPOLs
============================================================

WaterworldPipeline integra todos los módulos del pipeline a través de mixins:

    DiskSetupMixin         (disk_setup.py)
      → setup_grid()
      → setup_star()
      → initialize_simulation()


    DiskChemistryMixin     (disk_chemistry.py)
      → add_volatile_components()
      → add_ice_sigma_fields()

    SnowlinePhysicsMixin   (snowline_physics.py)
      → setup_physics()
      → add_snowline_fields()

    PressureBumpsMixin     (pressure_bumps.py)
      → setup_gap_kanagawa()
      → setup_gap_duffell()
      → reset_gap()

Uso típico
----------
    pipeline = WaterworldPipeline("mi_run/")
    pipeline.active_species = ["H2O", "CO2", "CO"]
    pipeline.setup_grid(rmin=1*c.au, rmax=300*c.au, Nr=200)
    pipeline.setup_star(M_star_Msun=1.0)
    pipeline.initialize_simulation()
    pipeline.add_volatile_components()
    pipeline.setup_physics()
    pipeline.add_snowline_fields()
    pipeline.add_ice_sigma_fields()
    pipeline.sim.update()
    pipeline.run_integration(t_end_years=1e6, num_snapshots=50)
"""

import math
import numpy as np
import dustpy.constants as c
from tripodpy import Simulation

from pipeline_methods import (
    DiskSetupMixin,
    DiskChemistryMixin,
    SnowlinePhysicsMixin,
    PressureBumpsMixin,
)


class WaterworldPipeline(
    DiskSetupMixin,
    DiskChemistryMixin,
    SnowlinePhysicsMixin,
    PressureBumpsMixin,
):
    """
    Pipeline maestro para simular la evolución del disco protoplanetario,
    snowlines (H2O, CO2, CO) y acumulación de pebbles usando tripodpy.

    Todos los métodos de configuración están en módulos mixin:
      - disk_setup.py       → grilla, estrella, inicialización, evolución estelar
      - disk_chemistry.py   → componentes volátiles, silicatos, fields de Σ
      - snowline_physics.py → v_frag(T), posiciones de snowline

    Parámetros configurables
    ------------------------
    active_species : list[str]
        Especies a inyectar como componentes (deben existir en chem.txt).
        Default: ["H2O", "CO2", "CO"].

    vfrag_params : dict
        {especie: (T_sub [K], v_frag [cm/s])} para el perfil v_frag(T).

    vfrag_silicates : float
        v_frag baseline de silicatos [cm/s]. Default: 100.

    M_star_Msun, R_star_Rsun, T_star_K : float
        Parámetros estelares en unidades solares / Kelvin.
    """

    def __init__(
        self, 
        datadir="output_pipeline",
        active_species=None,
        grid_rmin=0.5 * c.au,
        grid_rmax=100.0 * c.au,
        Nr=300,
        M_star_Msun=1.0,
        R_star_Rsun=2.1,  # ~2.1 R_sun para T Tauri con L=1 L_sun y T=4000K
        T_star_K=4000.0,   # T Tauri star
        gap_positions_au=None,
        snowline_au=2.7,
        alpha_gas=1e-3,
        M_disk_Msun=0.05,  # en masas solares, luego se convierte cgs
        v_frag_m_s=10.0,   # velocidad de fragmentación [m/s]
    ):
        self.sim     = Simulation()
        self.datadir = datadir
        
        self.alpha_gas = alpha_gas
        self.M_disk_Msun = M_disk_Msun
        self.v_frag_m_s = v_frag_m_s

        # ── Especies activas ──────────────────────────────────────────────────
        self.active_species = active_species if active_species is not None else ["H2O", "CO2", "CO"]

        # ── Parámetros estelares ──────────────────────────────────────────────
        self.M_star_Msun = M_star_Msun      # [M_sun]
        self.R_star_Rsun = R_star_Rsun      # [R_sun]  
        self.T_star_K    = T_star_K         # [K]      

        # ── Velocidades de fragmentación del hielo ────────────────────────────
        self.vfrag_params = {
            "H2O": (150.0, v_frag_m_s * 100.0),   # (T_sub [K], v_frag [cm/s])
            "CO2": ( 70.0,  500.0),
            "CO" : ( 25.0,  300.0),
        }
        self.vfrag_silicates = 100.0   # baseline refractario [cm/s]

        # ── Parámetros de gap planetario (PressureBumpsMixin) ─────────────────
        self.gap_alpha_ref   = 1e-3               # α de referencia [adim]
        self.gap_a_planet_au = 5.0                # semieje del planeta [AU]
        self.gap_M_planet    = 30. * c.M_earth    # masa del planeta [g — CGS]

        # ── Configuración automática de la simulación ─────────────────────────
        # 1. Configuración de malla y estrella
        self.setup_grid(rmin=grid_rmin, rmax=grid_rmax, Nr=Nr)
        self.setup_star(M_star_Msun=self.M_star_Msun, R_star_Rsun=self.R_star_Rsun, T_star_K=self.T_star_K)

        # 2. Refinamiento automático de la grilla (si el usuario provee posiciones)
        if gap_positions_au is not None:
            self.setup_refined_grid(gap_positions_au=gap_positions_au, snowline_au=snowline_au)

        # 3. Inicialización de los parámetros del gas base (previo a initialize)
        self.sim.ini.gas.alpha = self.alpha_gas
        self.sim.ini.gas.Mdisk = self.M_disk_Msun * c.M_sun

        # 4. Inicialización del motor base
        self.initialize_simulation()

        # 5. Química y Física de snowlines automáticas
        self.add_volatile_components()
        self.add_snowline_fields()
        self.setup_physics()
        self.add_ice_sigma_fields()
        
        # 6. Sincronización final pre-run
        self.sim.update()


    # ══════════════════════════════════════════════════════════════════════════
    # Integración temporal
    # ══════════════════════════════════════════════════════════════════════════

    def run_integration(self, t_end_years=1e5, num_snapshots=30, t_start_years=1e3):
        """
        Ejecuta la integración temporal y guarda snapshots en HDF5.

        Parameters
        ----------
        t_end_years : float
            Tiempo final de integración [yr]. Default: 1e5.
        num_snapshots : int
            Número de snapshots distribuidos log-uniformemente. Default: 30.
        t_start_years : float
            Primer snapshot [yr]. Default: 1e3.
        """
        print(f"Configurando integración:  "
              f"{t_start_years:.1e} → {t_end_years:.1e} yr  |  "
              f"{num_snapshots} snapshots")

        if t_start_years <= 0:
            raise ValueError(
                f"t_start_years debe ser > 0 (se usa log10). "
                f"Usa t_start_years=1e3 para empezar en 1000 yr."
            )
        self.sim.t.snapshots = np.logspace(
            math.log10(t_start_years),
            math.log10(t_end_years),
            num=num_snapshots,
        ) * c.year

        self.sim.writer.datadir   = self.datadir
        self.sim.writer.overwrite = True

        print("Iniciando integración...")
        self.sim.run()


# ══════════════════════════════════════════════════════════════════════════════
# Punto de entrada — ejemplo de uso completo
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    # Inicialización automática de todo el disco y física
    pipeline = WaterworldPipeline(
        datadir="data/post_pipeline/pipeline_v4",
        active_species=["H2O", "CO2", "CO"],
        grid_rmin=1 * c.au,
        grid_rmax=300 * c.au,
        Nr=200,
        M_star_Msun=1.0,
        R_star_Rsun=2.08,
        T_star_K=4000.0,
    )

    # Si quisiéramos agregar un gap, lo haríamos aquí:
    # pipeline.setup_gap_duffell(M_planet=317.8*c.M_earth, a_planet_au=5.0)
    # pipeline.sim.update()

    # Ejecutar integración
    pipeline.run_integration(t_end_years=1e6, num_snapshots=50)

    print("\nPipeline completado.")

"""
pipeline_snowlines.py — Pipeline principal del modelo PPOLs
============================================================

WaterworldPipeline integra todos los módulos del pipeline a través de mixins:

    DiskSetupMixin         (disk_setup.py)
      → setup_grid()
      → setup_star()
      → initialize_simulation()
      → setup_star_evolution()

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
    pipeline.setup_star_evolution()
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

    def __init__(self, datadir="output_pipeline"):
        self.sim     = Simulation()
        self.datadir = datadir

        # ── Parámetros estelares ──────────────────────────────────────────────
        self.M_star_Msun = 1.0      # [M_sun]
        self.R_star_Rsun = 2.0      # [R_sun]  radio joven inflado (Baraffe+ 2015)
        self.T_star_K    = 5778.0   # [K]      temperatura efectiva solar

        # ── Especies activas ──────────────────────────────────────────────────
        # Modificar ANTES de add_volatile_components().
        # Deben estar definidas en chem.txt.
        #   ["H2O"]               → solo agua
        #   ["H2O", "CO2"]        → agua + CO2
        #   ["H2O", "CO2", "CO"]  → completo (default)
        self.active_species = ["H2O", "CO2", "CO"]

        # ── Velocidades de fragmentación del hielo ────────────────────────────
        # Modificar ANTES de setup_physics() para cambiar la física de colisiones.
        # Referencias:
        #   Silicatos (baseline) : Birnstiel et al. (2012) A&A 539, A148
        #   H2O ice              : Gundlach & Blum (2015)  ApJ 798, 34
        #   CO2 ice              : Gundlach et al. (2018)  MNRAS 479, 1273
        #   CO  ice              : Dominik & Tielens (1997) ApJ 480, 647
        self.vfrag_params = {
            "H2O": (150.0, 1000.0),   # (T_sub [K], v_frag [cm/s])
            "CO2": ( 70.0,  500.0),
            "CO" : ( 25.0,  300.0),
        }
        self.vfrag_silicates = 100.0   # baseline refractario [cm/s]

        # ── Parámetros de gap planetario (PressureBumpsMixin) ─────────────────
        # Modificar ANTES de setup_gap_kanagawa() o setup_gap_duffell().
        # O pasar directamente como argumentos a cada método.
        #
        # IMPORTANTE: gap_alpha_ref es un escalar fijo capturado en el closure
        # del updater. NO se usa sim.gas.alpha para evitar circularidad.
        #
        # M_planet se da en CGS (gramos). Ejemplos:
        #   1.0 * c.M_jup    → Júpiter (~gap profundo)
        #   0.3 * c.M_jup    → Saturno (~gap intermedio)
        #   30. * c.M_earth  → ~0.1 M_Jup (pressure bump parcial)  ← default
        #   2.  * c.M_earth  → Super-Tierra (bump muy suave)
        self.gap_alpha_ref   = 1e-3               # α de referencia [adim]
        self.gap_a_planet_au = 5.0                # semieje del planeta [AU]
        self.gap_M_planet    = 30. * c.M_earth    # masa del planeta [g — CGS]


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
    pipeline = WaterworldPipeline("data/post_pipeline/pipeline_v4")

    # 1. Grilla y estrella
    pipeline.active_species = ["H2O", "CO2", "CO"]
    pipeline.setup_grid(rmin=1 * c.au, rmax=300 * c.au, Nr=200)
    pipeline.setup_star()                  # 1 M☉, 2 R☉, 5778 K

    # 2. Inicializar tripodpy
    pipeline.initialize_simulation()

    # 3. Química del disco
    pipeline.add_volatile_components()     # H2O, CO2, CO + silicatos

    # 4. Física de snowlines
    pipeline.setup_physics()               # v_frag(T) multi-especie
    pipeline.setup_star_evolution()        # contracción pre-MS → snowlines migran
    pipeline.add_snowline_fields()         # rsnow_H2O/CO2/CO → HDF5

    # 5. Fields de composición en HDF5
    pipeline.add_ice_sigma_fields()        # SigmaDust/SigmaGas por componente

    # 6. Sincronización y run
    pipeline.sim.update()
    pipeline.run_integration(t_end_years=1e6, num_snapshots=50)

    print("\nPipeline completado.")

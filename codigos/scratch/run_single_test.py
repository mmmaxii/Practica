import os
import sys
sys.path.append(os.path.abspath('.'))
import time
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

BASE_DIR = "data/5myr_test"

def run_test():
    label = "single_sup_earth_5au_no_star_evol"
    datadir = os.path.join(BASE_DIR, label)

    print(f"\n{'='*70}")
    print(f"[RUNNING TEST] {label}")
    print(f"[DATADIR] {datadir}")
    print(f"{'='*70}")
    t0 = time.time()

    pipeline = WaterworldPipeline(datadir)
    pipeline.active_species = ["H2O"]

    pipeline.setup_grid(rmin=1*c.au, rmax=100*c.au, Nr=200)
    pipeline.setup_star()

    # Refinamiento de grilla
    pipeline.setup_refined_grid(
        gap_positions_au=[5.0],
        num_gap=3,
        num_snow=2,
    )

    pipeline.initialize_simulation()
    pipeline.add_volatile_components()
    pipeline.setup_physics()

    # DESACTIVADO: Evolución estelar para ver migración natural
    # pipeline.setup_star_evolution()

    pipeline.add_snowline_fields()
    pipeline.add_ice_sigma_fields()

    # Estructura del disco
    pipeline.setup_gap_duffell(
        M_planet=5.0 * c.M_earth,
        a_planet_au=5.0,
    )

    pipeline.sim.update()

    # Integración
    pipeline.run_integration(
        t_start_years=100,
        t_end_years=5e6,
        num_snapshots=60,
    )

    elapsed = (time.time() - t0) / 60
    print(f"✓ '{label}' terminada en {elapsed:.1f} min.")

if __name__ == "__main__":
    os.makedirs(BASE_DIR, exist_ok=True)
    run_test()

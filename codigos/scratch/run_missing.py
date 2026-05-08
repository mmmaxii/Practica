"""
run_missing.py — Corre solo las simulaciones sinusoidales que faltan.
Detectadas comparando data_1e5/ con el plan de run_batch_1e5.py.
"""
import os
import time
import traceback
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

BASE_DIR = "data_1e5"

MISSING_RUNS = []

# ── Sinusoidal media que faltan ────────────────────────────────────────────────
for n_bumps in [30]:
    MISSING_RUNS.append({
        "label": f"sinusoidal_inner_amp_media_{n_bumps}gaps",
        "type": "sinusoidal",
        "amplitude": 3.0,
        "n_bumps": n_bumps,
        "r_inner_au": 1.0,
        "r_outer_au": 20.0
    })

for n_bumps in [20, 30]:
    MISSING_RUNS.append({
        "label": f"sinusoidal_outer_amp_media_{n_bumps}gaps",
        "type": "sinusoidal",
        "amplitude": 3.0,
        "n_bumps": n_bumps,
        "r_inner_au": 5.0,
        "r_outer_au": 50.0
    })

# ── Sinusoidal fuerte (amplitude=5.0) — todas faltan ──────────────────────────
extents = [
    {"name": "inner", "rin": 1.0,  "rout": 20.0},
    {"name": "outer", "rin": 5.0,  "rout": 50.0},
]
for ext in extents:
    for n_bumps in [5, 10, 15, 20, 30]:
        MISSING_RUNS.append({
            "label": f"sinusoidal_{ext['name']}_amp_fuerte_{n_bumps}gaps",
            "type": "sinusoidal",
            "amplitude": 5.0,
            "n_bumps": n_bumps,
            "r_inner_au": ext["rin"],
            "r_outer_au": ext["rout"]
        })


def run_simulation(cfg):
    label   = cfg["label"]
    datadir = os.path.join(BASE_DIR, label)

    # Skip si ya existe y tiene datos
    if os.path.isdir(datadir) and any(f.endswith(".hdf5") for f in os.listdir(datadir)):
        print(f"[SKIP] {label} ya existe con datos.")
        return

    print(f"\n{'='*70}\n[RUNNING] {label}\n[DATADIR] {datadir}\n{'='*70}")
    t0 = time.time()

    pipeline = WaterworldPipeline(datadir)
    pipeline.active_species = ["H2O"]
    pipeline.setup_grid(rmin=1*c.au, rmax=100*c.au, Nr=100)
    pipeline.setup_star()
    pipeline.initialize_simulation()
    pipeline.add_volatile_components()
    pipeline.setup_physics()
    pipeline.add_snowline_fields()
    pipeline.add_ice_sigma_fields()

    pipeline.setup_alpha_sinusoidal(
        amplitude   = cfg["amplitude"],
        n_bumps     = cfg["n_bumps"],
        r_inner_au  = cfg["r_inner_au"],
        r_outer_au  = cfg["r_outer_au"]
    )

    pipeline.sim.update()
    pipeline.run_integration(
        t_start_years = 100.0,
        t_end_years   = 1e5,
        num_snapshots = 50
    )
    print(f"✓ '{label}' terminada en {(time.time()-t0)/60:.1f} min.")


if __name__ == "__main__":
    os.makedirs(BASE_DIR, exist_ok=True)
    print(f"Simulaciones pendientes: {len(MISSING_RUNS)}")
    for i, cfg in enumerate(MISSING_RUNS, 1):
        print(f"\n--- [{i}/{len(MISSING_RUNS)}] ---")
        try:
            run_simulation(cfg)
        except Exception:
            print(f"[ERROR] '{cfg['label']}' falló:")
            traceback.print_exc()
            print("Saltando...\n")

    print(f"\n{'='*70}\n¡LOTE PENDIENTE TERMINADO!\nDatos en: {BASE_DIR}/")

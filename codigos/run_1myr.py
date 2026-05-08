"""
run_1myr.py — Batch de simulaciones a 1 Myr (1e6 años) con snowline dinámico.
==============================================================================

Configuraciones incluidas:
  - Gaps individuales: 5 tipos de planeta × 6 posiciones (1,2,3,5,7,10 AU)
  - Multi-planetas: 4 configuraciones representativas
  - Sinusoidal:
      3 inner (amp=0.5,1.0,2.0, n_bumps=10, r=[1,20] AU)
      3 outer  (amp=0.5,1.0,2.0, n_bumps=10, r=[5,50] AU)

Diferencias vs run_batch_1e5.py:
  - t_end = 1e6 yr  (10× más largo)
  - setup_star_evolution() activado → snowline migra durante el Myr
  - Snapshots: 60 (mejor resolución temporal)
  - Salida en data/1myr/
"""

import os
import time
import traceback
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

BASE_DIR = "data/1myr"

# ══════════════════════════════════════════════════════════════════════════════
# Definición de runs
# ══════════════════════════════════════════════════════════════════════════════

RUNS = []

# ── 1. Gaps individuales ──────────────────────────────────────────────────────
masses = {
    "sup_jup":   3.0 * 317.8 * c.M_earth,
    "jup":       317.8 * c.M_earth,
    "sat":       95.16 * c.M_earth,
    "nep":       17.15 * c.M_earth,
    "sup_earth": 5.0 * c.M_earth,
}
positions_1myr = [5.0, 7.0, 10.0, 15.0, 20.0]  # AU

for m_name, m_val in masses.items():
    for pos in positions_1myr:
        RUNS.append({
            "label":      f"single_{m_name}_{pos}au",
            "type":       "duffell",
            "M_planet":   m_val,
            "a_planet_au": pos,
        })

# ── 2. Multi-planetas (representativos) ──────────────────────────────────────
RUNS.extend([
    {
        "label": "multi_jup5_hjup7_jup10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["jup"],        "a_planet_au": 5.0},
            {"M_planet": 0.5 * masses["jup"],  "a_planet_au": 7.0},
            {"M_planet": masses["jup"],        "a_planet_au": 10.0},
        ]
    },
    {
        "label": "multi_sat3_jup5_sat7",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["sat"], "a_planet_au": 3.0},
            {"M_planet": masses["jup"], "a_planet_au": 5.0},
            {"M_planet": masses["sat"], "a_planet_au": 7.0},
        ]
    },
    {
        "label": "multi_nep3_sat5_jup10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["nep"], "a_planet_au": 3.0},
            {"M_planet": masses["sat"], "a_planet_au": 5.0},
            {"M_planet": masses["jup"], "a_planet_au": 10.0},
        ]
    },
    {
        "label": "multi_4jup_3_5_10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["jup"], "a_planet_au": 3.0},
            {"M_planet": masses["jup"], "a_planet_au": 5.0},
            {"M_planet": masses["jup"], "a_planet_au": 10.0},
        ]
    },
])

# ── 3. Sinusoidal — amplitudes bajas/medias ───────────────────────────────────
# inner: r = [1, 20] AU
# outer: r = [5, 50] AU
# amplitudes: 0.5 (muy suave), 1.0 (suave), 2.0 (media)
# n_bumps fijo = 10 (compromiso resolución / complejidad)
amplitudes_1myr = [0.5, 1.0, 2.0]
extents_1myr = [
    {"name": "middle", "rin": 5.0,  "rout": 20.0},
]

for ext in extents_1myr:
    for amp in amplitudes_1myr:
        amp_label = f"{amp:.1f}".replace(".", "p")   # "0p5", "1p0", "2p0"
        RUNS.append({
            "label":      f"sinusoidal_{ext['name']}_amp{amp_label}_10gaps",
            "type":       "sinusoidal",
            "amplitude":  amp,
            "n_bumps":    10,
            "r_inner_au": ext["rin"],
            "r_outer_au": ext["rout"],
        })

# ── Resumen ───────────────────────────────────────────────────────────────────
print(f"\nTotal de configuraciones: {len(RUNS)}")
print(f"  - Gaps individuales:  {5 * len(positions_1myr)}")
print(f"  - Multi-planetas:     4")
print(f"  - Sinusoidal:         {len(amplitudes_1myr) * len(extents_1myr)}")


# ══════════════════════════════════════════════════════════════════════════════
# Lógica de ejecución
# ══════════════════════════════════════════════════════════════════════════════

def run_simulation(cfg):
    label   = cfg["label"]
    datadir = os.path.join(BASE_DIR, label)

    # Skip si ya hay datos completos (50+ archivos HDF5)
    if os.path.isdir(datadir):
        hdfs = [f for f in os.listdir(datadir) if f.endswith(".hdf5")]
        if len(hdfs) >= 50:
            print(f"[SKIP] {label} ya tiene {len(hdfs)} snapshots.")
            return

    print(f"\n{'='*70}")
    print(f"[RUNNING] {label}")
    print(f"[DATADIR] {datadir}")
    print(f"{'='*70}")
    t0 = time.time()

    # ── Inicializar pipeline ───────────────────────────────────────────────
    pipeline = WaterworldPipeline(datadir)
    pipeline.active_species = ["H2O"]

    pipeline.setup_grid(rmin=1*c.au, rmax=100*c.au, Nr=200)
    pipeline.setup_star()

    # ── Refinamiento de grilla ANTES de initialize (patrón dustpylib) ────
    if cfg["type"] == "duffell":
        gap_pos = [cfg["a_planet_au"]]
    elif cfg["type"] == "duffell_multi":
        gap_pos = [p["a_planet_au"] for p in cfg["planets"]]
    else:  # sinusoidal — sin posición fija de gap
        gap_pos = None

    pipeline.setup_refined_grid(
        gap_positions_au = gap_pos,
        num_gap          = 3,   # niveles de refinamiento en el gap
        num_snow         = 2,   # niveles de refinamiento en la snowline H2O
    )

    pipeline.initialize_simulation()

    pipeline.add_volatile_components()
    pipeline.setup_physics()

    # DESACTIVADO: Snowline dinámico — para ver migración térmica pura
    # pipeline.setup_star_evolution()

    pipeline.add_snowline_fields()
    pipeline.add_ice_sigma_fields()

    # ── Instalar estructura del disco ──────────────────────────────────────
    if cfg["type"] == "duffell":
        pipeline.setup_gap_duffell(
            M_planet     = cfg["M_planet"],
            a_planet_au  = cfg["a_planet_au"],
        )
    elif cfg["type"] == "duffell_multi":
        pipeline.setup_gap_duffell_multi(planets=cfg["planets"])
    elif cfg["type"] == "sinusoidal":
        pipeline.setup_alpha_sinusoidal(
            amplitude   = cfg["amplitude"],
            n_bumps     = cfg["n_bumps"],
            r_inner_au  = cfg["r_inner_au"],
            r_outer_au  = cfg["r_outer_au"],
        )

    pipeline.sim.update()

    # ── Integración 1e3 → 1e6 yr con 60 snapshots (log-uniforme) ──────────
    pipeline.run_integration(
        t_start_years = 100,
        t_end_years   = 1e6,
        num_snapshots = 60,
    )


    elapsed = (time.time() - t0) / 60
    print(f"✓ '{label}' terminada en {elapsed:.1f} min.")


# ══════════════════════════════════════════════════════════════════════════════
# Main
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    os.makedirs(BASE_DIR, exist_ok=True)
    print(f"\nSalida en: {BASE_DIR}/")
    print(f"Total: {len(RUNS)} simulaciones a 1 Myr\n")

    for i, cfg in enumerate(RUNS, 1):
        print(f"\n--- Progreso: [{i}/{len(RUNS)}] ---")
        try:
            run_simulation(cfg)
        except Exception:
            print(f"\n[ERROR] '{cfg['label']}' falló:")
            traceback.print_exc()
            print("Saltando...\n")

    print(f"\n{'='*70}")
    print(f"¡LOTE 1 MYR COMPLETO!")
    print(f"Datos en: {BASE_DIR}/")

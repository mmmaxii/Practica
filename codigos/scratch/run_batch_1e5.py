import os
import time
import traceback
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

BASE_DIR = "data/1e5yr"

# ══════════════════════════════════════════════════════════════════════════════
# Configuración del Batch (aprox ~30 simulaciones)
# ══════════════════════════════════════════════════════════════════════════════

RUNS = []
# RUNS = [
#     {"label": "baseline", "type": "none"}
# ]

# 1. Gaps individuales (Super-Júpiter, Júpiter, Saturno, Neptuno, Super-Tierra)
masses = {
    "sup_jup":   3.0 * 317.8 * c.M_earth,
    "jup":       317.8 * c.M_earth,
    "sat":       95.16 * c.M_earth,
    "nep":       17.15 * c.M_earth,
    "sup_earth": 5.0 * c.M_earth
}
# Más densidad y extensión de posiciones (9 posiciones)
positions = [1.0, 1.5, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0]

for m_name, m_val in masses.items():
    for pos in positions:
        RUNS.append({
            "label": f"single_{m_name}_{pos}au",
            "type": "duffell",
            "M_planet": m_val,
            "a_planet_au": pos
        })

# 2. Combinaciones múltiples de gaps 
# (Expandido para ver casos de enjambres, gigantes múltiples, y cadencias)
RUNS.extend([
    {
        "label": "multi_jup5_hjup7_jup10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["jup"], "a_planet_au": 5.0},
            {"M_planet": 0.5 * masses["jup"], "a_planet_au": 7.0},
            {"M_planet": masses["jup"], "a_planet_au": 10.0}
        ]
    },
    {
        "label": "multi_sat3_jup5_sat7",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["sat"], "a_planet_au": 3.0},
            {"M_planet": masses["jup"], "a_planet_au": 5.0},
            {"M_planet": masses["sat"], "a_planet_au": 7.0}
        ]
    },
    {
        "label": "multi_nep3_sat5_jup10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["nep"], "a_planet_au": 3.0},
            {"M_planet": masses["sat"], "a_planet_au": 5.0},
            {"M_planet": masses["jup"], "a_planet_au": 10.0}
        ]
    },
    {
        "label": "multi_4jup_2_5_10_15",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["jup"], "a_planet_au": 2.0},
            {"M_planet": masses["jup"], "a_planet_au": 5.0},
            {"M_planet": masses["jup"], "a_planet_au": 10.0},
            {"M_planet": masses["jup"], "a_planet_au": 15.0}
        ]
    },
    {
        "label": "multi_mix_jup_sat_nep_3_6_9",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["jup"], "a_planet_au": 3.0},
            {"M_planet": masses["sat"], "a_planet_au": 6.0},
            {"M_planet": masses["nep"], "a_planet_au": 9.0}
        ]
    },
    {
        "label": "multi_swarm_5super_earths_2_to_10",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["sup_earth"], "a_planet_au": 2.0},
            {"M_planet": masses["sup_earth"], "a_planet_au": 4.0},
            {"M_planet": masses["sup_earth"], "a_planet_au": 6.0},
            {"M_planet": masses["sup_earth"], "a_planet_au": 8.0},
            {"M_planet": masses["sup_earth"], "a_planet_au": 10.0}
        ]
    },
    {
        "label": "multi_3sup_jup_5_10_20",
        "type": "duffell_multi",
        "planets": [
            {"M_planet": masses["sup_jup"], "a_planet_au": 5.0},
            {"M_planet": masses["sup_jup"], "a_planet_au": 10.0},
            {"M_planet": masses["sup_jup"], "a_planet_au": 20.0}
        ]
    }
])

# 3. Gaps sinusoidales ampliados 
amplitudes = {"suave": 1.0, "media": 3.0, "fuerte": 5.0, "extrema": 10.0}
n_bumps_list = [5, 10, 15, 20, 30]
extents = [
    {"name": "inner", "rin": 1.0, "rout": 20.0},
    {"name": "outer", "rin": 5.0, "rout": 50.0}
]

for amp_name, amp_val in amplitudes.items():
    for n_bumps in n_bumps_list:
        for ext in extents:
            RUNS.append({
                "label": f"sinusoidal_{ext['name']}_amp_{amp_name}_{n_bumps}gaps",
                "type": "sinusoidal",
                "amplitude": amp_val,
                "n_bumps": n_bumps,
                "r_inner_au": ext["rin"],
                "r_outer_au": ext["rout"]
            })


# ══════════════════════════════════════════════════════════════════════════════
# Lógica de Ejecución
# ══════════════════════════════════════════════════════════════════════════════

def run_simulation(cfg):
    label = cfg["label"]
    datadir = os.path.join(BASE_DIR, label)
    
    print(f"\n{'='*70}\n[RUNNING] {label}\n[DATADIR] {datadir}\n{'='*70}")
    t0 = time.time()
    
    pipeline = WaterworldPipeline(datadir)
    
    # Química con todos los hielos
    pipeline.active_species = ["H2O"]
    
    # Setups dictados
    pipeline.setup_grid(rmin=1*c.au, rmax=100*c.au, Nr=100)
    pipeline.setup_star()
    pipeline.initialize_simulation()
    
    pipeline.add_volatile_components()
    pipeline.setup_physics()
    
    # SIN evolución de la estrella (Desactivado).
    # pipeline.setup_star_evolution() 
    
    pipeline.add_snowline_fields()
    pipeline.add_ice_sigma_fields()
    
    # Instalar gap correspondiente (ESTÁTICO A T=0 SOLO duffellos)
    if cfg["type"] == "duffell":
        pipeline.setup_gap_duffell(M_planet=cfg["M_planet"], a_planet_au=cfg["a_planet_au"])
    elif cfg["type"] == "duffell_multi":
        pipeline.setup_gap_duffell_multi(planets=cfg["planets"])
    elif cfg["type"] == "sinusoidal":
        pipeline.setup_alpha_sinusoidal(
            amplitude=cfg["amplitude"],
            n_bumps=cfg["n_bumps"],
            r_inner_au=cfg["r_inner_au"],
            r_outer_au=cfg["r_outer_au"]
        )
        
    pipeline.sim.update()
    
    # 1e3 a 1e5 (100,000 años), 50 snapshots
    pipeline.run_integration(
        t_start_years = 100.0,
        t_end_years   = 1e5,
        num_snapshots = 50
    )
    
    print(f"✓ '{label}' terminada exitosamente en {(time.time()-t0)/60:.1f} minutos.")


if __name__ == "__main__":
    os.makedirs(BASE_DIR, exist_ok=True)
    print(f"Total de configuraciones a procesar: {len(RUNS)}")
    
    for i, cfg in enumerate(RUNS, 1):
        print(f"\n--- Progreso General: [{i}/{len(RUNS)}] ---")
        try:
            run_simulation(cfg)
        except Exception as e:
            print(f"\n[ERROR CRÍTICO] La run '{cfg['label']}' falló:")
            traceback.print_exc()
            print("Saltando a la siguiente configuración...\n")

    print(f"\n{'='*70}\n¡LOTE COMPLETO TERMINADO!\nTodos los datos están en: {BASE_DIR}/")

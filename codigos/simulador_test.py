import sys
import os
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

if len(sys.argv) < 2:
    print("Error: Se requiere indicar el tipo de prueba ('normal' o 'gap')")
    sys.exit(1)

tipo_prueba = sys.argv[1]
os.makedirs(os.path.join("runs_geryon", "test"), exist_ok=True)

if tipo_prueba == "normal":
    print("=== Corriendo Prueba Normal (Sin Gap) ===")
    ruta = os.path.join("runs_geryon", "test", "run_normal")
    pipeline = WaterworldPipeline(
        datadir=ruta,
        active_species=[],
        grid_rmin=0.5 * c.au,
        grid_rmax=100.0 * c.au,
        Nr=300,
        M_star_Msun=1.0,
        R_star_Rsun=2.1,
        T_star_K=4000.0,
        alpha_gas=1e-3,
        M_disk_Msun=0.05
    )
    pipeline.run_integration(t_start_years=100, t_end_years=10000, num_snapshots=10)

elif tipo_prueba == "gap":
    print("=== Corriendo Prueba con Gap (10 au, 1 M_jup) ===")
    ruta = os.path.join("runs_geryon", "test", "run_gap_10au_1Mjup")
    pipeline = WaterworldPipeline(
        datadir=ruta,
        active_species=[],
        grid_rmin=0.5 * c.au,
        grid_rmax=100.0 * c.au,
        Nr=300,
        M_star_Msun=1.0,
        R_star_Rsun=2.1,
        T_star_K=4000.0,
        gap_positions_au=[10.0],
        alpha_gas=1e-3,
        M_disk_Msun=0.05
    )
    pipeline.setup_gap_duffell(M_planet=1.0 * c.M_jup, a_planet_au=10.0)
    pipeline.sim.update()
    pipeline.run_integration(t_start_years=100, t_end_years=10000, num_snapshots=10)

else:
    print(f"Tipo de prueba desconocido: {tipo_prueba}")
    sys.exit(1)
    
print("Prueba finalizada.")

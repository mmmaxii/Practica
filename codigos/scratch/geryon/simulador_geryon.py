import sys
import os
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

# 1. ATRAPAR LOS 3 PARÁMETROS DESDE LA CONSOLA
if len(sys.argv) < 4:
    print("Error: Se requieren 3 parámetros: r_gap, M_gap, alpha")
    sys.exit(1)

radio_gap = float(sys.argv[1])
masa_planeta = float(sys.argv[2])
viscosidad = float(sys.argv[3])

job_name = f"run_r{radio_gap}_m{masa_planeta}_a{viscosidad}"
ruta_guardado = os.path.join("runs_geryon", "10Myr", job_name)

print("==================================================")
print(f"Iniciando configuración de disco estructurado:")
print(f"-> R_gap : {radio_gap} AU")
print(f"-> M_gap : {masa_planeta} M_jup")
print(f"-> Alpha : {viscosidad}")
print("==================================================")

# ====================================================
# CONFIGURACIÓN DEL PIPELINE
# ====================================================
pipeline = WaterworldPipeline(
    datadir=ruta_guardado,
    active_species=[],
    grid_rmin=0.5 * c.au,
    grid_rmax=100.0 * c.au,
    Nr=300,
    M_star_Msun=1.0,
    R_star_Rsun=2.1,
    T_star_K=4000.0,
    gap_positions_au=[radio_gap], # Refina grilla en el gap
    alpha_gas=viscosidad,
    M_disk_Msun=0.05
)

# Configurar el gap con el método DUFFELL
pipeline.setup_gap_duffell(M_planet=masa_planeta * c.M_jup, a_planet_au=radio_gap)
pipeline.sim.update()

# Correr de 100 yr a 10 Myr con 100 snapshots
pipeline.run_integration(t_start_years=100, t_end_years=1e7, num_snapshots=100)

print(f"Simulación finalizada. Datos guardados en {ruta_guardado}")

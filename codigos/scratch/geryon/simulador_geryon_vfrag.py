import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

# 1. ATRAPAR LOS 4 PARÁMETROS DESDE LA CONSOLA
if len(sys.argv) < 5:
    print("Error: Se requieren 4 parámetros: r_gap, M_gap, alpha, v_frag")
    sys.exit(1)

radio_gap = float(sys.argv[1])
masa_planeta = float(sys.argv[2])
viscosidad = float(sys.argv[3])
vel_frag = float(sys.argv[4])

job_name = f"run_r{radio_gap}_m{masa_planeta}_a{viscosidad}_v{vel_frag}"
ruta_guardado = os.path.join("runs_geryon", f"v_frag_{int(vel_frag)}ms", job_name)

print("==================================================")
print(f"Iniciando configuración de disco estructurado:")
print(f"-> R_gap  : {radio_gap} AU")
print(f"-> M_gap  : {masa_planeta} M_jup")
print(f"-> Alpha  : {viscosidad}")
print(f"-> V_frag : {vel_frag} m/s")
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
    M_disk_Msun=0.05,
    v_frag_m_s=vel_frag
)

# Configurar el gap con el método DUFFELL
pipeline.setup_gap_duffell(M_planet=masa_planeta * c.M_jup, a_planet_au=radio_gap)
pipeline.sim.update()

# Correr de 100 yr a 10 Myr con 100 snapshots
pipeline.run_integration(t_start_years=100, t_end_years=1e7, num_snapshots=100)

print(f"Simulación finalizada. Datos guardados en {ruta_guardado}")

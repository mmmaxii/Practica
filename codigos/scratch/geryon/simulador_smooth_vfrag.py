# -*- coding: utf-8 -*-
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

if len(sys.argv) < 3:
    print("Error: Se requieren 2 parámetros: alpha, v_frag")
    sys.exit(1)

viscosidad = float(sys.argv[1])
vel_frag = float(sys.argv[2])

job_name = f"run_smooth_a{viscosidad}_v{int(vel_frag)}"
ruta_guardado = os.path.join("runs_geryon", f"vf_{int(vel_frag)}ms", "smooth_baselines", job_name)

print("==================================================")
print(f"Iniciando configuración de disco SMOOTH:")
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
    alpha_gas=viscosidad,
    M_disk_Msun=0.05,
    v_frag_m_s=vel_frag
)

# NO llamamos a setup_gap_duffell, por lo que el disco no tiene subestructuras.
pipeline.sim.update()

# Correr de 100 yr a 10 Myr con 100 snapshots
pipeline.run_integration(t_start_years=100, t_end_years=1e7, num_snapshots=100)

print(f"Simulación finalizada. Datos guardados en {ruta_guardado}")

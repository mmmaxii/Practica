# -*- coding: utf-8 -*-
import sys
import os

# Asegurar que Python encuentre el módulo base en la raíz del proyecto
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

if len(sys.argv) < 7:
    print("Error: Se requieren 6 parámetros: n_gaps, amplitud, alpha, v_frag, nombre_modelo, r_min")
    sys.exit(1)

n_bumps = int(sys.argv[1])
amplitud = float(sys.argv[2])
alpha_val = float(sys.argv[3])
v_frag_val = float(sys.argv[4])
nombre_modelo = sys.argv[5]
r_min_val = float(sys.argv[6])

job_name = f"run_ngap{n_bumps}_A{amplitud}_a{alpha_val}_rmin{r_min_val}"
ruta_guardado = os.path.join("runs_geryon", nombre_modelo, job_name)

print("==================================================")
print(f"Iniciando simulación Sinusoidal Grid:")
print(f"N_gaps  : {n_bumps}")
print(f"Amplitud: {amplitud}")
print(f"Alpha   : {alpha_val}")
print(f"V_frag  : {v_frag_val} m/s")
print(f"R_inner : 10.0 AU (Fijo)")
print(f"R_outer : 100.0 AU (Fijo)")
print(f"R_min   : {r_min_val} AU (Dinámico)")
print(f"Destino : {ruta_guardado}")
print("==================================================")

os.makedirs(ruta_guardado, exist_ok=True)

pipeline = WaterworldPipeline(
    datadir=ruta_guardado,
    active_species=[], 
    grid_rmin=r_min_val * c.au, 
    grid_rmax=100.0 * c.au,
    Nr=200,
    M_star_Msun=1.0,
    R_star_Rsun=2.1,
    T_star_K=4000.0,
    alpha_gas=alpha_val,
    M_disk_Msun=0.05,
    v_frag_m_s=v_frag_val
)

# Imprimir el perfil de alpha sinusoidal en ventana inmutable [10, 100] AU
pipeline.setup_alpha_sinusoidal(
    alpha_ref=alpha_val,
    amplitude=amplitud,
    n_bumps=n_bumps,
    r_inner_au=10.0,
    r_outer_au=100.0,
    imprint=True
)
pipeline.sim.update()

# Integracion principal (10 Myr, 100 snapshots)
pipeline.run_integration(t_start_years=100, t_end_years=1e7, num_snapshots=100)

print(f"Simulación {job_name} finalizada.")

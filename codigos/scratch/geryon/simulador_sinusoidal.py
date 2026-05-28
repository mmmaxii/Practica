# -*- coding: utf-8 -*-
import sys
import os
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

if len(sys.argv) < 6:
    print("Error: Se requieren 5 parámetros: n_bumps, amplitud, r_inner, r_outer, nombre_modelo")
    sys.exit(1)

n_bumps = int(sys.argv[1])
amplitud = float(sys.argv[2])
r_inner = float(sys.argv[3])
r_outer = float(sys.argv[4])
nombre_modelo = sys.argv[5]

viscosidad = 1e-3

job_name = f"{nombre_modelo}_A{amplitud}"
ruta_guardado = os.path.join("runs_geryon", "sinusoidales", job_name)

print("==================================================")
print(f"Iniciando simulación Sinusoidal - Modelo {nombre_modelo}:")
print(f"N_bumps : {n_bumps}")
print(f"Amplitud: {amplitud}")
print(f"R_inner : {r_inner} AU")
print(f"R_outer : {r_outer} AU")
print("==================================================")

pipeline = WaterworldPipeline(
    datadir=ruta_guardado,
    active_species=[], 
    grid_rmin=0.7 * c.au, 
    grid_rmax=100.0 * c.au,
    Nr=200,
    M_star_Msun=1.0,
    R_star_Rsun=2.1,
    T_star_K=4000.0,
    alpha_gas=viscosidad,
    M_disk_Msun=0.05
)

# Imprimir el perfil de alpha sinusoidal
pipeline.setup_alpha_sinusoidal(
    alpha_ref=viscosidad,
    amplitude=amplitud,
    n_bumps=n_bumps,
    r_inner_au=r_inner,
    r_outer_au=r_outer,
    imprint=True
)
pipeline.sim.update()

# Integracion principal
pipeline.run_integration(t_start_years=100, t_end_years=1e7, num_snapshots=100)

print(f"Simulación {job_name} finalizada.")

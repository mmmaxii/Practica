# -*- coding: utf-8 -*-
import sys
import os
import numpy as np
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

# 1. Parámetros
radio_gap = 15.0
masa_planeta = 1.0 # M_jup
viscosidad = 1e-3
tiempo_insercion_gap = 1e5 # 100.000 años
tiempo_final = 1e6         # 1 millon de años

job_name = f"run_retrasado_r{radio_gap}_m{masa_planeta}_tgap{tiempo_insercion_gap}"
ruta_guardado = os.path.join("runs_geryon", job_name)

print("==================================================")
print(f"Iniciando simulación con GAP RETRASADO:")
print(f"R_gap : {radio_gap} AU")
print(f"M_gap : {masa_planeta} M_jup")
print(f"Alpha : {viscosidad}")
print(f"Tiempo de inserción: {tiempo_insercion_gap} años")
print("==================================================")

# ====================================================
# CONFIGURACIÓN DEL PIPELINE
# ====================================================
pipeline = WaterworldPipeline(
    datadir=ruta_guardado,
    active_species=[], # Agrega ["H2O", "silicates"] si necesitas la química
    grid_rmin=1 * c.au,
    grid_rmax=100.0 * c.au,
    Nr=100,
    M_star_Msun=1.0,
    R_star_Rsun=2.1,
    T_star_K=4000.0,
    gap_positions_au=[radio_gap], # Refina grilla en el gap
    alpha_gas=viscosidad,
    M_disk_Msun=0.05
)

# NOTA: NO insertamos el gap todavía. El disco empieza "liso" (smooth).

# ====================================================
# FASE 1: Evolución sin Gap (100 yr -> tiempo_insercion_gap)
# ====================================================
print(f"\n--- FASE 1: Evolución de disco liso hasta {tiempo_insercion_gap} años ---")
snapshots_fase1 = 30
tiempos_fase1 = np.logspace(np.log10(100), np.log10(tiempo_insercion_gap), snapshots_fase1)

pipeline.sim.t.snapshots = tiempos_fase1 * c.year
pipeline.sim.writer.datadir = ruta_guardado
pipeline.sim.writer.overwrite = True 

pipeline.sim.run()

# ====================================================
# FASE 2: Insertar el Gap Planetario
# ====================================================
print(f"\n--- FASE 2: Insertando el planeta de {masa_planeta} M_jup ---")
pipeline.setup_gap_duffell(M_planet=masa_planeta * c.M_jup, a_planet_au=radio_gap)

# Imprimir explícitamente el gap en las densidades actuales para que el hueco 
# se cree instantáneamente en el gas y el polvo:
pipeline._imprint_gap(pipeline.sim.gas.alpha, pipeline.sim.gas.alpha.copy(), label="gap retrasado") 
pipeline.sim.update()

# ====================================================
# FASE 3: Evolución con Gap (tiempo_insercion_gap -> 10 Myr)
# ====================================================
print(f"\n--- FASE 3: Evolución con gap hasta {tiempo_final} años ---")
snapshots_fase2 = 20
# Tomamos desde el tiempo del gap hasta el final (evitando repetir el primer snapshot)
tiempos_fase2 = np.logspace(np.log10(tiempo_insercion_gap), np.log10(tiempo_final), snapshots_fase2)[1:]

# DustPy siempre empieza a numerar desde data0000.hdf5 al llamar a sim.run()
# Para no sobreescribir la Fase 1, guardamos en una carpeta temporal y luego renombramos.
ruta_fase2_temp = os.path.join(ruta_guardado, "fase2_temp")
pipeline.sim.writer.datadir = ruta_fase2_temp
pipeline.sim.writer.overwrite = True

pipeline.sim.t.snapshots = tiempos_fase2 * c.year
pipeline.sim.run()

# Unificar los archivos
import shutil
print("\nReordenando y unificando snapshots HDF5...")
archivos_fase2 = sorted([f for f in os.listdir(ruta_fase2_temp) if f.startswith("data") and f.endswith(".hdf5")])

indice_actual = snapshots_fase1  # El siguiente número disponible después de la Fase 1
for archivo in archivos_fase2:
    nombre_viejo = os.path.join(ruta_fase2_temp, archivo)
    nombre_nuevo = os.path.join(ruta_guardado, f"data{indice_actual:04d}.hdf5")
    shutil.move(nombre_viejo, nombre_nuevo)
    indice_actual += 1

# Limpiar carpeta temporal
os.rmdir(ruta_fase2_temp)

print(f"\nSimulación completada. Datos unificados y guardados en {ruta_guardado}")

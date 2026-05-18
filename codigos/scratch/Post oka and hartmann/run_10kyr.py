import numpy as np
import dustpy.constants as c
import shutil
import os

from pipeline_snowlines import WaterworldPipeline

output_dir = "output_10kyr_oka"

# Limpiar directorio si existe
if os.path.exists(output_dir):
    shutil.rmtree(output_dir)

print("Inicializando pipeline con física de snowline (Oka 2011)...")
pipeline = WaterworldPipeline(
    active_species=[],  # Array vacío sin componentes extra
    M_star_Msun=1.0,
    R_star_Rsun=2.1,
    T_star_K=4000.0,
    datadir=output_dir
)

# Ejecutar simulación de 100 a 10,000 años con 10 snapshots
pipeline.run_integration(t_end_years=1000000.0, num_snapshots=80, t_start_years=100.0)

print("\n---------------------------------------------------------")
print("Evolución del Snowline - Resultados Extraídos del HDF5:")
print("---------------------------------------------------------")

from dustpy import hdf5writer

# Leer los datos guardados por la simulación
wrtr = hdf5writer()
wrtr.datadir = output_dir
data = wrtr.read.all()

times_yr = data.t[:] / c.year
r_snow_au = data.dust.r_snow[:] / c.au

print("Tiempo [yr]     |   r_snow [AU]")
print("--------------------------------")
for t, r in zip(times_yr, r_snow_au):
    print(f"{t:14.2f}  |  {r:10.2f}")

# -*- coding: utf-8 -*-
import os
import subprocess

# Caso de prueba exigente (prometedor)
N = 7
A = 1.0
alpha = 0.001
v_frag = 10.0
nombre_modelo = "benchmarks"

# Valores de r_min a testear
r_mins = [0.5, 0.7]

print("Generando scripts PBS para el Benchmark A/B de r_min...")

for r_min in r_mins:
    v_str = str(int(v_frag)) if v_frag.is_integer() else str(v_frag)
    
    # Nombre de directorio y job
    job_name = f"run_ngap{N}_A{A}_a{alpha}_rmin{r_min}"
    
    # Nombre corto para PBS (max 15 chars)
    pbs_job_name = f"BM_rmin{r_min}"
    
    pbs_content = f"""#PBS -S /bin/bash
#PBS -V
#PBS -N {pbs_job_name}
#PBS -j eo
#PBS -M m.valderrama.geryon2@gmail.com
#PBS -m a
#PBS -l select=1:ppn=1
#PBS -l walltime=48:00:00

cd $PBS_O_WORKDIR

# Cargar compilador y entorno
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

# Ejecutar simulador (N_gaps, Amplitud, Alpha, V_frag, Nombre_Modelo, r_min)
python simulador_sinusoidal_grid.py {N} {A} {alpha} {v_frag} {nombre_modelo} {r_min}
"""
    
    script_file = f"submit_tmp_{pbs_job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Descomentar si se corre en entorno local vs geryon. En Geryon, ejecutar subprocess.
    subprocess.run(["qsub", script_file])
    os.remove(script_file)
    
    print(f"Enviado Benchmark: {job_name} con r_min = {r_min}")

print("Terminado. Los 2 jobs de benchmark enviados a Geryon.")

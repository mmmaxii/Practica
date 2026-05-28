# -*- coding: utf-8 -*-
import os
import itertools
import subprocess

# 1. CONFIGURACIÓN DE GAPS PLANETARIOS
r_gap = [5, 10, 15, 20]
M_gap = [0.01, 0.1, 0.3, 0.5, 1, 3, 5]  # Masas de Júpiter
alphas = [0.005]

# Crear carpeta de resultados principal
os.makedirs(os.path.join("runs_geryon", "10Myr"), exist_ok=True)

# 2. BUCLE COMBINATORIO
n_jobs = len(r_gap) * len(M_gap) * len(alphas)
for r, m, a in itertools.product(r_gap, M_gap, alphas):
    
    # Nombre identificador (ej: run_r10_m1_a0.001)
    job_name = f"run_r{r}_m{m}_a{a}"
    
    pbs_content = f"""#PBS -S /bin/bash
#PBS -V
#PBS -N {job_name}
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

# Ejecutar el simulador enviando los 3 parámetros en orden: radio, masa, alpha
python simulador_geryon.py {r} {m} {a}
"""
    
    script_file = f"submit_{job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Enviar al clúster
    subprocess.run(["qsub", script_file])
    os.remove(script_file)

print(f"Las {n_jobs} configuraciones de gaps han sido enviadas a la cola para alpha {alphas[0]}")

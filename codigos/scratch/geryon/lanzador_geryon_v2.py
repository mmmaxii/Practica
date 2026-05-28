# -*- coding: utf-8 -*-
import os
import itertools
import subprocess

# 1. CONFIGURACIÓN DE GAPS PLANETARIOS (Focalizado en baja masa)
r_gap = [5, 10, 12, 15, 17, 20, 25, 30]
M_gap = [0.01, 0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3]  # Masas de Júpiter
alphas = [1e-3]

# Crear carpeta de resultados principal
os.makedirs(os.path.join("runs_geryon", "10Myr_around0.1"), exist_ok=True)

# 2. BUCLE COMBINATORIO (Generará 72 jobs)
for r, m, a in itertools.product(r_gap, M_gap, alphas):
    
    # Nombre identificador (ej: run_r10_m0.05_a0.001)
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

# Ejecutar el simulador v2 enviando los 3 parámetros en orden: radio, masa, alpha
python simulador_geryon_v2.py {r} {m} {a}
"""
    
    script_file = f"submit_{job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Enviar al clúster
    subprocess.run(["qsub", script_file])
    os.remove(script_file)

print("Las 72 configuraciones de gaps (M_gap <= 0.3 M_jup) han sido enviadas a la cola")

# -*- coding: utf-8 -*-
import os
import itertools
import subprocess

# 1. CONFIGURACIÓN DE GAPS PLANETARIOS Y V_FRAG
r_gap = [5, 10, 15, 20]
M_gap = [0.01, 0.05, 0.1, 0.15, 0.3, 0.5, 1, 2, 3]  # Masas de Júpiter
alphas = [0.0001, 0.0003, 0.0005, 0.0007]
v_frags = [1]

# 2. BUCLE COMBINATORIO
for r, m, a, v in itertools.product(r_gap, M_gap, alphas, v_frags):
    
    # Crear carpeta de resultados principal si no existe
    os.makedirs(os.path.join("runs_geryon", f"v_frag_{int(v)}ms"), exist_ok=True)

    # Nombre identificador
    job_name = f"run_r{r}_m{m}_a{a}_v{v}"
    
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

# Ejecutar el simulador enviando los 4 parámetros en orden: radio, masa, alpha, v_frag
python simulador_geryon_vfrag.py {r} {m} {a} {v}
"""
    
    script_file = f"submit_{job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Enviar al clúster
    subprocess.run(["qsub", script_file])
    os.remove(script_file)

total_jobs = len(r_gap) * len(M_gap) * len(alphas) * len(v_frags)
print(f"Las {total_jobs} configuraciones para v_frag = 1m/s han sido enviadas a la cola.")

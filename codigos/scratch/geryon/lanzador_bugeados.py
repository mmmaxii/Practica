# -*- coding: utf-8 -*-
import os
import subprocess

# Lista explícita de configuraciones que se quedaron pegadas (r_gap, M_gap, alpha)
configs = [
    (5, 3.0, 0.001),
    (10, 0.01, 0.001),
    (15, 0.01, 0.001),
    (15, 5.0, 0.001),
    (20, 0.01, 0.001),
]

for r, m, a in configs:
    # Nombre identificador original
    # Ej: run_r10_m0.01_a0.001
    job_name = f"run_r{r}_m{m}_a{a}"
    
    pbs_content = f"""#PBS -S /bin/bash
#PBS -V
#PBS -N {job_name}
#PBS -j eo
#PBS -M m.valderrama.geryon2@gmail.com
#PBS -m a
#PBS -l select=1:ppn=1
#PBS -l walltime=72:00:00

cd $PBS_O_WORKDIR

# Cargar compilador y entorno
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

# Ejecutar el simulador normal apuntando a la carpeta de 10Myr original
python simulador_geryon.py {r} {m} {a}
"""
    
    script_file = f"submit_{job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Enviar al clúster
    subprocess.run(["qsub", script_file])
    os.remove(script_file)

print(f"Se han re-enviado los {len(configs)} jobs pendientes con 72 hrs de walltime.")

# -*- coding: utf-8 -*-
import os
import subprocess

# Configuración (sin gaps)
alphas = [1e-4, 1e-3, 1e-2]

# Crear carpeta de resultados
os.makedirs(os.path.join("runs_geryon", "10Myr_smooth"), exist_ok=True)

# Bucle para lanzar los 3 jobs (uno por cada alpha)
for a in alphas:
    job_name = f"run_smooth_a{a}"
    
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

# Ejecutar el simulador smooth
python simulador_smooth.py {a}
"""
    
    script_file = f"submit_{job_name}.pbs"
    with open(script_file, "w") as file:
        file.write(pbs_content)
    
    # Enviar al clúster
    subprocess.run(["qsub", script_file])
    os.remove(script_file)

print(f"Se han enviado {len(alphas)} configuraciones de discos smooth a la cola")

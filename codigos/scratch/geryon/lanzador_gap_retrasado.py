# -*- coding: utf-8 -*-
import os
import subprocess

masas = [0.01, 0.1, 0.3, 0.5, 1, 3, 5]
tiempos_insercion = [1e5, 5e5, 1e6, 2e6, 5e6]

os.makedirs(os.path.join("runs_geryon", "delayed_gaps"), exist_ok=True)
total_jobs = len(masas) * len(tiempos_insercion)

print(f"Generando scripts PBS para {total_jobs} configuraciones de gaps retrasados...")

for m in masas:
    for t in tiempos_insercion:
        job_name = f"dgap_m{m}_t{t:.0f}"
        
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

# Ejecutar simulador: masa_planeta, tiempo_insercion_gap
python simulador_gap_retrasado.py {m} {t:.0f}
"""
        
        script_file = f"submit_{job_name}.pbs"
        with open(script_file, "w") as file:
            file.write(pbs_content)
        
        subprocess.run(["qsub", script_file])
        os.remove(script_file)
        
        print(f"Enviado a la cola: {script_file}")

print(f"Terminado. {total_jobs} jobs enviados a Geryon.")

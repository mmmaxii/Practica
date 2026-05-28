# -*- coding: utf-8 -*-
import os
import subprocess

# Lista maestra de amplitudes a explorar
amplitudes_explorar = [0.5, 0.7, 1.0, 2.0, 3.0, 5.0]

# CONFIGURACIÓN DE LOS CASOS SINUSOIDALES
casos = [
    {
        "nombre": "MuchosGaps", 
        "n_bumps": 7, 
        "r_inner": 10.0, 
        "r_outer": 100.0, 
        "amplitudes": amplitudes_explorar
    },
    {
        "nombre": "Intermedio", 
        "n_bumps": 5, 
        "r_inner": 15.0, 
        "r_outer": 100.0, 
        "amplitudes": amplitudes_explorar
    },
    {
        "nombre": "PocosGaps", 
        "n_bumps": 3, 
        "r_inner": 20.0, 
        "r_outer": 100.0, 
        "amplitudes": amplitudes_explorar
    }
]

os.makedirs(os.path.join("runs_geryon", "sinusoidales"), exist_ok=True)
total_jobs = sum([len(c["amplitudes"]) for c in casos])

print(f"Generando scripts PBS para {total_jobs} configuraciones...")

for caso in casos:
    n_bumps = caso["n_bumps"]
    r_inner = caso["r_inner"]
    r_outer = caso["r_outer"]
    
    for A in caso["amplitudes"]:
        job_name = f"sinus_{caso['nombre']}_A{A}"
        
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

# Ejecutar simulador
python simulador_sinusoidal.py {n_bumps} {A} {r_inner} {r_outer} {caso['nombre']}
"""
        
        script_file = f"submit_{job_name}.pbs"
        with open(script_file, "w") as file:
            file.write(pbs_content)
        
        subprocess.run(["qsub", script_file])
        os.remove(script_file)
        
        print(f"Enviado a la cola: {script_file}")

print(f"Terminado. {total_jobs} jobs enviados a Geryon.")

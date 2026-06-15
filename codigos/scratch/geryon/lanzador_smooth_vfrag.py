# -*- coding: utf-8 -*-
import os
import subprocess

# Configuración (sin gaps)
v_frags = [1, 3, 10]
alphas = [1e-4, 3e-4, 5e-4, 7e-4, 1e-3, 3e-3, 5e-3, 1e-2]

total_jobs = len(v_frags) * len(alphas)

print("==================================================")
print(f"Preparando lanzador para discos SMOOTH (sin gap)")
print(f"Total de configuraciones a simular: {total_jobs}")
print("==================================================")

count = 0
for v in v_frags:
    # Crear carpeta de resultados para este v_frag
    os.makedirs(os.path.join("runs_geryon", f"vf_{int(v)}ms", "smooth_baselines"), exist_ok=True)
    
    for a in alphas:
        count += 1
        job_name = f"run_smooth_a{a}_v{int(v)}"
        
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

# Ejecutar el simulador smooth con v_frag
python simulador_smooth_vfrag.py {a} {v}
"""
        
        script_file = f"submit_{job_name}.pbs"
        with open(script_file, "w") as file:
            file.write(pbs_content)
        
        # Enviar al clúster (solo funciona en Geryon)
        try:
            subprocess.run(["qsub", script_file], check=True)
            print(f"[{count}/{total_jobs}] Job enviado: {job_name}")
        except FileNotFoundError:
            # Si no estamos en el clúster (e.g. corriendo en Windows local)
            print(f"[{count}/{total_jobs}] PBS generado: {script_file} (no se pudo ejecutar qsub)")
        
        # Eliminamos el archivo .pbs local para no ensuciar
        if os.path.exists(script_file):
            os.remove(script_file)

print("\n¡Proceso de generación finalizado!")

# -*- coding: utf-8 -*-
import os
import subprocess

# Parámetros de la grilla
N_gaps = [3, 5, 7, 10, 15, 20]
amplitudes = [0.5, 0.7, 1.0, 2.0, 3.0, 5.0]
alphas = [0.0005, 0.001]
v_frags = [10.0, 3.0, 1.0]

total_jobs = len(N_gaps) * len(amplitudes) * len(alphas) * len(v_frags)
print(f"Generando scripts PBS para {total_jobs} configuraciones...")

for v in v_frags:
    v_str = str(int(v)) if v.is_integer() else str(v)
    for a in alphas:
        for A in amplitudes:
            for n in N_gaps:
                nombre_modelo = f"sinusoidal/vf_{v_str}"
                job_name = f"run_ngap{n}_A{A}_a{a}"
                
                # Nombre corto para PBS (max 15 chars permitidos por PBS)
                pbs_job_name = f"S{n}A{A}a{a}v{v_str}"
                
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
python simulador_sinusoidal_grid.py {n} {A} {a} {v} {nombre_modelo} 0.7
"""
                
                script_file = f"submit_tmp_{pbs_job_name}.pbs"
                with open(script_file, "w") as file:
                    file.write(pbs_content)
                
                # Descomentar para enviar a Geryon directamente
                subprocess.run(["qsub", script_file])
                os.remove(script_file)
                
                print(f"Enviado: {job_name} para vfrag {v_str}")

print(f"Terminado. {total_jobs} jobs enviados a Geryon.")

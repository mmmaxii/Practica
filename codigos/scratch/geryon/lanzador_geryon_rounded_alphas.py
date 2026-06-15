# -*- coding: utf-8 -*-
import os
import subprocess

# 1. CONFIGURACIÓN DE GAPS PLANETARIOS ESPECÍFICA POR ALPHA ('rounded')
r_gap_list = [5, 10, 12, 15, 17, 20, 25, 30]

# Definimos el rango del Sweet Spot que le corresponde a cada turbulencia
m_gap_por_alpha = {
    0.0001: [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1],      # Muy bajo alpha -> Gaps pequeños
    0.0005: [0.03, 0.05, 0.08, 0.1, 0.15, 0.2, 0.25, 0.3],        # Similar al 0.001
    0.003:  [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]         # Alto alpha -> Sweet Spot desplazado a gaps masivos
}

v_frag = 10 # v_frag estricto a 10 m/s

# Puntos de la grilla general que ya tenemos calculados
r_gen_completos = [5, 10, 15, 20]
m_gen_completos = [0.01, 0.05, 0.1, 0.15, 0.3, 0.5, 1.0, 2.0, 3.0]

# 2. BUCLE COMBINATORIO
total_jobs = 0
for alpha, m_gaps in m_gap_por_alpha.items():
    for r in r_gap_list:
        for m in m_gaps:
            
            # Evitar repetir cálculos que ya existen en tu carpeta vf_10ms/general
            if r in r_gen_completos and m in m_gen_completos:
                continue
            
            # Nombre identificador
            job_name = f"run_rounded_r{r}_m{m}_a{alpha}_v{v_frag}"
            
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

# Ejecutar el simulador enviando 5 parámetros: radio, masa, alpha, v_frag, carpeta_destino
python simulador_geryon_vfrag.py {r} {m} {alpha} {v_frag} rounded_a{alpha}
"""
            
            script_file = f"submit_{job_name}.pbs"
            with open(script_file, "w") as file:
                file.write(pbs_content)
            
            # Enviar al clúster
            subprocess.run(["qsub", script_file])
            os.remove(script_file)
            total_jobs += 1

print(f"Las {total_jobs} configuraciones de gaps 'rounded' han sido enviadas a la cola de Geryon con 48 hrs de límite.")

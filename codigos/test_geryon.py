import os
import subprocess

print("==================================================")
print("Iniciando lanzador de pruebas (Test Runs): 100 yr hasta 10000 yr")
print("==================================================")

os.makedirs(os.path.join("runs_geryon", "test"), exist_ok=True)

# 1. Prueba normal (sin gap)
job_name_1 = "test_normal"
pbs_content_1 = f"""#PBS -S /bin/bash
#PBS -V
#PBS -N {job_name_1}
#PBS -j eo
#PBS -l select=1:ppn=1
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

python simulador_test.py normal
"""

script_file_1 = f"submit_{job_name_1}.pbs"
with open(script_file_1, "w") as file:
    file.write(pbs_content_1)
subprocess.run(["qsub", script_file_1])
os.remove(script_file_1)

# 2. Prueba con gap (10 au, 1 M_jup)
job_name_2 = "test_gap"
pbs_content_2 = f"""#PBS -S /bin/bash
#PBS -V
#PBS -N {job_name_2}
#PBS -j eo
#PBS -l select=1:ppn=1
#PBS -l walltime=10:00:00

cd $PBS_O_WORKDIR
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

python simulador_test.py gap
"""

script_file_2 = f"submit_{job_name_2}.pbs"
with open(script_file_2, "w") as file:
    file.write(pbs_content_2)
subprocess.run(["qsub", script_file_2])
os.remove(script_file_2)

print("¡Las 2 pruebas han sido enviadas a la cola del clúster (Geryon)!")

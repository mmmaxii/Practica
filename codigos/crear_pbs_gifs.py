import os
import glob

# Rutas base
dest_base = "data/figures/gif"
os.makedirs(dest_base, exist_ok=True)

pbs_header = """#PBS -S /bin/bash
#PBS -V
#PBS -N gifs_nuevos
#PBS -j eo
#PBS -M m.valderrama.geryon2@gmail.com
#PBS -m abe
#PBS -l select=1:ppn=8
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

echo "Generando nuevos GIFs..."

"""

commands = []

# 1. 10Myr (solo nuevos alphas: 0.0005, 0.003, 0.005)
alphas_nuevos = ["0.0005", "0.003", "0.005"]
runs_10myr = glob.glob("data/runs/10Myr/run_*")
for rpath in runs_10myr:
    run_name = os.path.basename(rpath)
    # extraer alpha
    alpha_str = None
    for p in run_name.split('_'):
        if p.startswith('a'):
            alpha_str = p[1:]
            break
    if alpha_str in alphas_nuevos:
        dest_dir = os.path.join(dest_base, f"10Myr_a{alpha_str}")
        dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
        cmd = f"python disk_visualizer/disk_viz.py --datadir \"{rpath}\" --gif \"{dest_gif}\" --rmax 30 --every 2 --fps 15"
        commands.append((dest_dir, cmd))

# 2. 10Myr_round0.1
runs_round = glob.glob("data/runs/10Myr_round0.1/run_*")
for rpath in runs_round:
    run_name = os.path.basename(rpath)
    alpha_str = None
    for p in run_name.split('_'):
        if p.startswith('a'):
            alpha_str = p[1:]
            break
    if alpha_str:
        dest_dir = os.path.join(dest_base, f"10Myr_round0.1_a{alpha_str}")
        dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
        cmd = f"python disk_visualizer/disk_viz.py --datadir \"{rpath}\" --gif \"{dest_gif}\" --rmax 30 --every 2 --fps 15"
        commands.append((dest_dir, cmd))

# 3. Sinusoidal
runs_sinu = glob.glob("data/runs/10Myr Sinusoidal/*_A*")
for rpath in runs_sinu:
    if not os.path.isdir(rpath): continue
    run_name = os.path.basename(rpath)
    dest_dir = os.path.join(dest_base, "Sinusoidal")
    dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
    cmd = f"python disk_visualizer/disk_viz.py --datadir \"{rpath}\" --gif \"{dest_gif}\" --rmax 30 --every 2 --fps 15"
    commands.append((dest_dir, cmd))

# 4. Delayed gaps
runs_delay = glob.glob("data/runs/10Myr delayed/run_*")
for rpath in runs_delay:
    if not os.path.isdir(rpath): continue
    run_name = os.path.basename(rpath)
    dest_dir = os.path.join(dest_base, "Delayed")
    dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
    cmd = f"python disk_visualizer/disk_viz.py --datadir \"{rpath}\" --gif \"{dest_gif}\" --rmax 30 --every 2 --fps 15"
    commands.append((dest_dir, cmd))

# Crear carpetas de destino y armar archivo final
with open("submit_gifs_nuevos.pbs", "w", encoding="utf-8") as f:
    f.write(pbs_header)
    for dest_dir, cmd in commands:
        os.makedirs(dest_dir, exist_ok=True)
        # Reemplazar contrabarras de Windows por barras normales para bash
        cmd_linux = cmd.replace('\\', '/')
        f.write(cmd_linux + "\n")
        
print(f"Generado submit_gifs_nuevos.pbs con {len(commands)} comandos.")

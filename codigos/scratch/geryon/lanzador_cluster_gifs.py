import os
import glob
import subprocess
import sys

def collect_commands():
    dest_base = "data/figures/gif"
    os.makedirs(dest_base, exist_ok=True)
    commands = []

    # 1. 10Myr (solo nuevos alphas: 0.0005, 0.003, 0.005)
    alphas_nuevos = ["0.0005", "0.003", "0.005"]
    runs_10myr = glob.glob("runs_geryon/10Myr/run_*")
    for rpath in runs_10myr:
        if not os.path.isdir(rpath): continue
        run_name = os.path.basename(rpath)
        alpha_str = None
        for p in run_name.split('_'):
            if p.startswith('a'):
                alpha_str = p[1:]
                break
        if alpha_str in alphas_nuevos:
            dest_dir = os.path.join(dest_base, f"10Myr_a{alpha_str}")
            os.makedirs(dest_dir, exist_ok=True)
            dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
            rpath_linux = rpath.replace('\\', '/')
            dest_gif_linux = dest_gif.replace('\\', '/')
            commands.append(f"python disk_visualizer/disk_viz.py --datadir \"{rpath_linux}\" --gif \"{dest_gif_linux}\" --rmax 30 --every 2 --fps 15")

    # 2. 10Myr_around0.1
    runs_round = glob.glob("runs_geryon/10Myr_around0.1/run_*")
    for rpath in runs_round:
        if not os.path.isdir(rpath): continue
        run_name = os.path.basename(rpath)
        alpha_str = None
        for p in run_name.split('_'):
            if p.startswith('a'):
                alpha_str = p[1:]
                break
        if alpha_str:
            dest_dir = os.path.join(dest_base, f"10Myr_round0.1_a{alpha_str}")
            os.makedirs(dest_dir, exist_ok=True)
            dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
            rpath_linux = rpath.replace('\\', '/')
            dest_gif_linux = dest_gif.replace('\\', '/')
            commands.append(f"python disk_visualizer/disk_viz.py --datadir \"{rpath_linux}\" --gif \"{dest_gif_linux}\" --rmax 30 --every 2 --fps 15")

    # 3. Sinusoidales
    runs_sinu = glob.glob("runs_geryon/sinusoidales/*_A*")
    for rpath in runs_sinu:
        if not os.path.isdir(rpath): continue
        run_name = os.path.basename(rpath)
        dest_dir = os.path.join(dest_base, "Sinusoidal")
        os.makedirs(dest_dir, exist_ok=True)
        dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
        rpath_linux = rpath.replace('\\', '/')
        dest_gif_linux = dest_gif.replace('\\', '/')
        commands.append(f"python disk_visualizer/disk_viz.py --datadir \"{rpath_linux}\" --gif \"{dest_gif_linux}\" --rmax 30 --every 2 --fps 15")

    # 4. Delayed gaps
    runs_delay = glob.glob("runs_geryon/delayed_gaps/run_*")
    for rpath in runs_delay:
        if not os.path.isdir(rpath): continue
        run_name = os.path.basename(rpath)
        dest_dir = os.path.join(dest_base, "Delayed")
        os.makedirs(dest_dir, exist_ok=True)
        dest_gif = os.path.join(dest_dir, f"{run_name}.gif")
        rpath_linux = rpath.replace('\\', '/')
        dest_gif_linux = dest_gif.replace('\\', '/')
        commands.append(f"python disk_visualizer/disk_viz.py --datadir \"{rpath_linux}\" --gif \"{dest_gif_linux}\" --rmax 30 --every 2 --fps 15")

    return commands

def run_worker(cmd):
    # Ejecuta el bash command ocultando la salida estándar (para que el log de PBS no quede de miles de líneas)
    subprocess.run(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    # Parseo básico para imprimir el nombre del run que se completó
    try:
        gif_name = cmd.split("--gif ")[1].split(".gif")[0].split("/")[-1]
        print(f"Terminado: {gif_name}")
    except:
        pass

def main():
    # MODO WORKER: Esto es lo que ejecuta el PBS Job en el cluster
    if len(sys.argv) > 1 and sys.argv[1] == "--run-workers":
        import concurrent.futures
        commands = collect_commands()
        print(f"Iniciando ThreadPoolExecutor con 20 hilos para {len(commands)} tareas...")
        
        # ThreadPoolExecutor nativo de Python 3.11 lanzando 20 subprocesos en paralelo
        with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:
            executor.map(run_worker, commands)
            
        print("¡Todos los GIFs han sido generados exitosamente en paralelo!")
        return

    # MODO LANZADOR: Esto es lo que ejecutas tú en la consola
    commands = collect_commands()
    if not commands:
        print("No se encontraron corridas para generar GIFs. Revisa si las carpetas de datos existen.")
        return

    # PBS pidiendo 20 nucleos e incluyendo la notificacion por mail original
    pbs_header = """#PBS -S /bin/bash
#PBS -V
#PBS -N gifs_threads
#PBS -j eo
#PBS -M m.valderrama.geryon2@gmail.com
#PBS -m abe
#PBS -l select=1:ppn=20
#PBS -l walltime=24:00:00

cd $PBS_O_WORKDIR
module purge
module load devtoolset/devtoolset-11
source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311

echo "Ejecutando proceso maestro de multithreading con 20 nucleos..."
python lanzador_cluster_gifs.py --run-workers
"""
    pbs_filename = "submit_gifs_multithread.pbs"
    with open(pbs_filename, "w", encoding="utf-8", newline='\n') as f:
        f.write(pbs_header)
        
    print(f"Archivo {pbs_filename} creado exitosamente pidiendo 20 núcleos.")
    print(f"Se procesarán {len(commands)} GIFs en paralelo ocupando ThreadPoolExecutor de Python 3.11.")
    
    # Enviar al cluster
    print("Enviando el trabajo unificado al cluster con qsub...")
    try:
        subprocess.run(["qsub", pbs_filename], stdout=subprocess.DEVNULL, check=True)
        print("¡Trabajo enviado exitosamente a la cola!")
    except Exception as e:
        print(f"Error enviando el trabajo: {e}")

if __name__ == "__main__":
    main()

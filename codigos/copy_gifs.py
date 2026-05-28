import os
import shutil
import glob

source_dir = "data/runs/10Myr"
dest_base = "data/figures/gif"

runs = glob.glob(os.path.join(source_dir, "run_*"))

count = 0
for rpath in runs:
    run_name = os.path.basename(rpath)
    # Ejemplo run_name: run_r10.0_m0.01_a0.01
    parts = run_name.split('_')
    
    alpha_str = None
    for p in parts:
        if p.startswith('a'):
            alpha_str = p[1:]
            break
            
    if not alpha_str:
        continue
        
    gif_path = os.path.join(rpath, "gif", "disco.gif")
    if os.path.exists(gif_path):
        dest_dir = os.path.join(dest_base, f"10Myr_a{alpha_str}")
        os.makedirs(dest_dir, exist_ok=True)
        
        dest_path = os.path.join(dest_dir, f"{run_name}.gif")
        shutil.copy2(gif_path, dest_path)
        count += 1
        print(f"Copiado: {run_name} -> {dest_dir}")

print(f"Total de GIFs copiados: {count}")

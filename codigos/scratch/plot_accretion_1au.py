import os
import sys
sys.path.insert(0, os.path.abspath('.'))
import glob
# pyrefly: ignore [missing-import]
import numpy as np          
# pyrefly: ignore [missing-import]
import matplotlib.pyplot as plt
# pyrefly: ignore [missing-import]
import matplotlib.cm as cm
# pyrefly: ignore [missing-import]
import dustpy.constants as c
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

def plot_1au_mass_evolution(base_dir="data/1myr"):
    run_dirs = sorted(glob.glob(os.path.join(base_dir, "*")))
    run_dirs = [d for d in run_dirs if os.path.isdir(d)]
    
    if not run_dirs:
        print(f"No run directories found in {base_dir}")
        return

    plt.figure(figsize=(10, 7))
    
    colors = cm.viridis(np.linspace(0, 1, len(run_dirs)))
    
    for i, datadir in enumerate(run_dirs):
        label = os.path.basename(datadir)
        print(f"Processing {label}...")
        
        try:
            pam = PebbleAccretionModule3.from_datadir(datadir)
            results = pam.run_growth([1.0]) # Embrión en 1 AU
            
            hist = results[1.0]
            if len(hist) > 0:
                t_yr = hist[:, 0] / 3.156e7
                M_core = hist[:, 1] / c.M_earth
                
                plt.plot(t_yr, M_core, lw=1.5, label=label[:20], color=colors[i])
        except Exception as e:
            print(f"  Error processing {label}: {e}")

    plt.xlabel('Time [yr]', fontsize=12)
    plt.ylabel(r'Embryo Mass [$M_\oplus$]', fontsize=12)
    plt.title('Core Growth at 1 AU (1 Myr Runs)', fontsize=14)
    plt.grid(True, alpha=0.5)
    
    # Put legend outside
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=2)
    plt.tight_layout()
    
    out_file = os.path.join(base_dir, 'mass_evolution_1au.png')
    plt.savefig(out_file, dpi=300)
    plt.close()
    print(f"\nSaved plot to {out_file}")

if __name__ == "__main__":
    plot_1au_mass_evolution()

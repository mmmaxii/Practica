import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

# Constantes cgs
c_M_sun = 1.9891e33
c_year = 3.15576e7

base_dir = r"c:\astro\Codigos practica + docs + papers\codigos\data\runs\5myr"
runs_10au = [
    "single_nep_10.0au",
    "single_sat_10.0au",
    "single_jup_10.0au",
    "single_sup_jup_10.0au"
]
labels = [
    "Neptuno",
    "Saturno",
    "Júpiter",
    "Super Júpiter"
]

# Colores en gradiente
colors = cm.viridis(np.linspace(0.1, 0.9, len(runs_10au)))

plt.figure(figsize=(10, 6))

for run_idx, run_name in enumerate(runs_10au):
    datadir = os.path.join(base_dir, run_name)
    if not os.path.exists(datadir):
        print(f"Directory not found: {datadir}")
        continue
        
    snaps = [f for f in os.listdir(datadir) if f.startswith('data') and f.endswith('.hdf5')]
    snaps.sort()
    
    if not snaps:
        print(f"No hdf5 files in {datadir}")
        continue

    times_yr = []
    mdots_msun_yr = []
    
    for snap_file in snaps:
        filepath = os.path.join(datadir, snap_file)
        try:
            with h5py.File(filepath, 'r') as f:
                t_yr = f['t'][()] / c_year
                sigma = f['gas']['Sigma'][:]
                nu = f['gas']['nu'][:]
                
                mdot_array = 3 * np.pi * nu * sigma
                mdot_inner_cgs = mdot_array[0] 
                mdot_msun_yr = mdot_inner_cgs * (c_year / c_M_sun)
                
                times_yr.append(t_yr)
                mdots_msun_yr.append(mdot_msun_yr)
        except Exception as e:
            pass # ignore corrupted/incomplete snaps

    times_yr = np.array(times_yr)
    mdots_msun_yr = np.array(mdots_msun_yr)
    
    if len(times_yr) > 0:
        plt.plot(times_yr, mdots_msun_yr, '-', color=colors[run_idx], linewidth=2.5, label=f"{labels[run_idx]} (10 AU)")
        print(f"{labels[run_idx]} final Mdot: {mdots_msun_yr[-1]:.2e} M_sun/yr")

plt.xlabel("Tiempo [yr]", fontsize=12)
plt.ylabel(r"$\dot{M}_{\mathrm{gas}}$ [$M_\odot$/yr]", fontsize=12)
plt.title("Evolución temporal de la Tasa de Acreción de Gas\n(5 Myr, Planetas en 10 AU)", fontsize=14)
plt.grid(True, which="both", ls="--", alpha=0.7)
plt.yscale('log')
plt.xlim(left=max(1, times_yr[1]*0.5) if times_yr[0] == 0 else times_yr[0])
plt.legend(fontsize=11)
plt.tight_layout()

plot_path = "mdot_evolution_10au_gradient.png"
plt.savefig(plot_path, dpi=300)
print(f"Saved plot to {plot_path}")

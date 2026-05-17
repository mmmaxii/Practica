import h5py
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import os

# Constantes cgs
c_M_sun = 1.9891e33
c_year = 3.15576e7

base_dir = r"c:\astro\Codigos practica + docs + papers\codigos\data\runs\5myr"

gap_distances = ["1.0au", "3.0au", "5.0au", "7.0au", "10.0au"]
planet_types = ["single_sup_earth", "single_nep", "single_sat", "single_jup", "single_sup_jup"]
planet_labels = ["Super Tierra", "Neptuno", "Saturno", "Júpiter", "Super Júpiter"]

for dist in gap_distances:
    print(f"\n--- Procesando gap en {dist} ---")
    plt.figure(figsize=(10, 6))
    
    # We will only use colors for the planets that actually exist for this distance
    # Let's find out which exist first
    existing_planets = []
    for p_idx, p_type in enumerate(planet_types):
        datadir = os.path.join(base_dir, f"{p_type}_{dist}")
        if os.path.exists(datadir):
            existing_planets.append(p_idx)
            
    if not existing_planets:
        print(f"No hay simulaciones para {dist}.")
        plt.close()
        continue
        
    colors = cm.viridis(np.linspace(0.1, 0.9, len(existing_planets)))
    
    all_times = []
    
    for c_idx, p_idx in enumerate(existing_planets):
        p_type = planet_types[p_idx]
        p_label = planet_labels[p_idx]
        datadir = os.path.join(base_dir, f"{p_type}_{dist}")
        
        snaps = [f for f in os.listdir(datadir) if f.startswith('data') and f.endswith('.hdf5')]
        snaps.sort()
        
        if not snaps:
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
            plt.plot(times_yr, mdots_msun_yr, '-', color=colors[c_idx], linewidth=2.5, label=f"{p_label} ({dist})")
            print(f"{p_label} final Mdot: {mdots_msun_yr[-1]:.2e} M_sun/yr")
            if len(all_times) == 0:
                all_times = times_yr

    plt.xlabel("Tiempo [yr]", fontsize=12)
    plt.ylabel(r"$\dot{M}_{\mathrm{gas}}$ [$M_\odot$/yr]", fontsize=12)
    plt.title(f"Evolución temporal de la Tasa de Acreción de Gas\n(5 Myr, Planetas en {dist})", fontsize=14)
    plt.grid(True, which="both", ls="--", alpha=0.7)
    plt.yscale('log')
    
    if len(all_times) > 0:
        plt.xlim(left=max(1, all_times[1]*0.5) if all_times[0] == 0 else all_times[0])
        
    plt.legend(fontsize=11)
    plt.tight_layout()

    plot_path = f"mdot_evolution_{dist.replace('.0', '')}_gradient.png"
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved plot to {plot_path}")

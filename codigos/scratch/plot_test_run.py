import os
import sys
sys.path.insert(0, os.path.abspath('.'))
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import dustpy.constants as c
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

def plot_test_results(datadir="data/1myr_test/single_jup_5au_no_star_evol"):
    files = sorted(glob.glob(os.path.join(datadir, "data*.hdf5")))
    if not files:
        print(f"No HDF5 files found in {datadir}")
        return

    # 1. Temperature vs Distance (Snapshot 0)
    with h5py.File(files[0], 'r') as f:
        r_au = f['grid/r'][:] / c.au
        T_gas = f['gas/T'][:]
    
    plt.figure(figsize=(8, 6))
    plt.plot(r_au, T_gas, lw=2, color='firebrick')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Distance [AU]')
    plt.ylabel('Temperature [K]')
    plt.title('Gas Temperature vs Distance (Initial Snapshot)')
    plt.grid(True, which='both', ls='--', alpha=0.5)
    plt.axhline(150, color='blue', linestyle='--', label='H2O Snowline (150K)')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(datadir, 'temp_vs_distance_snap0.png'))
    plt.close()

    # 2. Snowline vs Time
    times_yr = []
    rsnow_au = []
    
    for fpath in files:
        with h5py.File(fpath, 'r') as f:
            times_yr.append(f['t'][()] / 3.156e7)
            if 'grid/rsnow_H2O' in f:
                rsnow_au.append(f['grid/rsnow_H2O'][()] / c.au)
            else:
                rsnow_au.append(np.nan)
                
    plt.figure(figsize=(8, 6))
    plt.plot(times_yr, rsnow_au, lw=2, color='blue')
    plt.xlabel('Time [yr]')
    plt.ylabel('H2O Snowline Position [AU]')
    plt.title('Snowline Migration Over Time')
    plt.grid(True, alpha=0.5)
    plt.tight_layout()
    plt.savefig(os.path.join(datadir, 'snowline_vs_time.png'))
    plt.close()

    # 3. Core Growth Over Time (using PA3Py)
    print("Running PebbleAccretion3 to get core growth...")
    try:
        pam = PebbleAccretionModule3.from_datadir(datadir)
        results = pam.run_growth([5.0]) # Embrión en 5 AU
        
        hist = results[5.0]
        if len(hist) > 0:
            t_pa = hist[:, 0] / 3.156e7
            M_core = hist[:, 1] / c.M_earth
            
            plt.figure(figsize=(8, 6))
            plt.plot(t_pa, M_core, lw=2, color='green')
            plt.xlabel('Time [yr]')
            plt.ylabel(r'Embryo Mass [$M_\oplus$]')
            plt.title('Core Growth Over Time (Embryo at 5 AU)')
            plt.grid(True, alpha=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(datadir, 'core_growth_vs_time.png'))
            plt.close()
            print(f"Plots saved in {datadir}")
    except Exception as e:
        print(f"Failed to run PA3Py: {e}")

if __name__ == "__main__":
    plot_test_results()

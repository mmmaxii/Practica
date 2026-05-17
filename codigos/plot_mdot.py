import simframe as sf
import matplotlib.pyplot as plt
import numpy as np
import os

# Constantes cgs
c_M_sun = 1.9891e33
c_year = 3.15576e7
c_au = 1.495978707e13

datadir = r"c:\astro\Codigos practica + docs + papers\codigos\data\gaps_pipeline\t_1e6_optimised\baseline"

# Count snapshots
snaps = [f for f in os.listdir(datadir) if f.startswith('data') and f.endswith('.hdf5')]
num_snaps = len(snaps)

times_yr = []
mdots_msun_yr = []

print(f"Reading {num_snaps} snapshots from {datadir}...")

for i in range(num_snaps):
    try:
        snap = sf.read_output(datadir, i)
        t_yr = snap.t / c_year
        r = snap.grid.r
        sigma = snap.gas.Sigma
        nu = snap.gas.nu
        
        mdot_array = 3 * np.pi * nu * sigma
        mdot_inner_cgs = mdot_array[0] # tomando el borde interno
        
        mdot_msun_yr = mdot_inner_cgs * (c_year / c_M_sun)
        
        times_yr.append(t_yr)
        mdots_msun_yr.append(mdot_msun_yr)
    except Exception as e:
        print(f"Error reading snapshot {i}: {e}")

times_yr = np.array(times_yr)
mdots_msun_yr = np.array(mdots_msun_yr)

plt.figure(figsize=(8, 5))
plt.plot(times_yr, mdots_msun_yr, 'o-', color='crimson')
plt.xlabel("Tiempo [yr]")
plt.ylabel(r"$\dot{M}_{\mathrm{gas}}$ [$M_\odot$/yr]")
plt.title("Tasa de Acreción de Gas en el Borde Interno (Baseline 1 Myr)")
plt.grid(True)
plt.yscale('log')
plt.tight_layout()
plt.savefig("mdot_evolution.png", dpi=300)
print(f"Final Mdot at {times_yr[-1]:.2e} yr: {mdots_msun_yr[-1]:.2e} M_sun/yr")
print("Saved plot to mdot_evolution.png")

# Also save the data to a text file for the user
data_to_save = np.column_stack((times_yr, mdots_msun_yr))
np.savetxt("mdot_evolution_data.txt", data_to_save, header="Time_yr Mdot_Msun_yr")
print("Saved data to mdot_evolution_data.txt")

import h5py
import matplotlib.pyplot as plt
import numpy as np
import os

# Constantes cgs
c_M_sun = 1.9891e33
c_year = 3.15576e7

datadir = r"c:\astro\Codigos practica + docs + papers\codigos\data\runs\5myr\single_nep_10.0au"

# Encontrar los archivos hdf5 y ordenarlos
snaps = [f for f in os.listdir(datadir) if f.startswith('data') and f.endswith('.hdf5')]
snaps.sort()
num_snaps = len(snaps)

times_yr = []
mdots_msun_yr = []

print(f"Reading {num_snaps} snapshots from {datadir}...")

for i, snap_file in enumerate(snaps):
    filepath = os.path.join(datadir, snap_file)
    try:
        with h5py.File(filepath, 'r') as f:
            t_yr = f['t'][()] / c_year
            sigma = f['gas']['Sigma'][:]
            nu = f['gas']['nu'][:]
            
            # M_dot = 3 * pi * nu * Sigma
            mdot_array = 3 * np.pi * nu * sigma
            
            # Tomamos el valor en el borde interno del disco
            mdot_inner_cgs = mdot_array[0] 
            
            mdot_msun_yr = mdot_inner_cgs * (c_year / c_M_sun)
            
            times_yr.append(t_yr)
            mdots_msun_yr.append(mdot_msun_yr)
    except Exception as e:
        print(f"Error reading {snap_file}: {e}")

times_yr = np.array(times_yr)
mdots_msun_yr = np.array(mdots_msun_yr)

# Graficar
plt.figure(figsize=(9, 6))
plt.plot(times_yr, mdots_msun_yr, 'o-', color='darkorange', linewidth=2, markersize=5)
plt.xlabel("Tiempo [yr]", fontsize=12)
plt.ylabel(r"$\dot{M}_{\mathrm{gas}}$ [$M_\odot$/yr]", fontsize=12)
plt.title("Evolución temporal de la Tasa de Acreción de Gas\n(5 Myr, Neptuno en 10 AU, Borde Interno)", fontsize=14)
plt.grid(True, which="both", ls="--", alpha=0.7)
plt.yscale('log')

# Limitar x axis si hay t=0 (log plot)
if len(times_yr) > 1 and times_yr[0] == 0:
    plt.xlim(left=max(1, times_yr[1]*0.5))

plt.tight_layout()
plot_path = "mdot_evolution_single_nep_10au.png"
plt.savefig(plot_path, dpi=300)

if len(times_yr) > 0:
    print(f"Final Mdot at {times_yr[-1]:.2e} yr: {mdots_msun_yr[-1]:.2e} M_sun/yr")
print(f"Saved plot to {plot_path}")

# Guardar los datos también
data_path = "mdot_evolution_data_single_nep_10au.txt"
data_to_save = np.column_stack((times_yr, mdots_msun_yr))
np.savetxt(data_path, data_to_save, header="Time_yr Mdot_Msun_yr", fmt="%.5e")
print(f"Saved data to {data_path}")

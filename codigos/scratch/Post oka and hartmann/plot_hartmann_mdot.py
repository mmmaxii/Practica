import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Rango temporal de 100 yr a 10^7 yr
t_yr = np.logspace(2, 7, 500)
t_myr = t_yr / 1e6

# Parámetros de Hartmann et al.
mdot_1myr = 1e-8 # M_sun / yr
etas = np.arange(1.5, 2.81, 0.1)

fig, ax_main = plt.subplots(figsize=(10, 7))

# Colormap para el gradiente de etas
colors = cm.viridis(np.linspace(0.1, 0.9, len(etas)))

for i, eta in enumerate(etas):
    mdot_t = mdot_1myr * (t_myr)**(-eta)
    # Solo agregar label para la leyenda si es el primero (1.5) o el último (2.8)
    if i == 0 or i == len(etas) - 1:
        line_label = fr'$\eta = {eta:.1f}$'
    else:
        line_label = '_nolegend_'
        
    ax_main.plot(np.log10(t_myr), np.log10(mdot_t), linewidth=2.0, color=colors[i], label=line_label)

# Líneas de referencia
ax_main.axvline(0.0, color='gray', linestyle='--', alpha=0.7, label='1 Myr')
ax_main.axhline(-8.0, color='gray', linestyle=':', alpha=0.7)

ax_main.set_xlabel(r"$\log_{10}(t \,[\mathrm{Myr}])$", fontsize=16)
ax_main.set_ylabel(r"$\log_{10}(\dot{M}_{\mathrm{gas}} \,[M_\odot/\mathrm{yr}])$", fontsize=16)
ax_main.set_title("T Tauri stars", fontsize=20, fontweight='bold')
ax_main.tick_params(axis='both', which='major', labelsize=14)
ax_main.grid(True, which="both", ls="--", alpha=0.7)

# Leyenda con solo las líneas etiquetadas (1.5, 2.8 y 1 Myr)
ax_main.legend(fontsize=14, loc='lower left')

# ==========================================
# Crear el Zoom (Inset Axes) arriba a la derecha
# ==========================================
axins = ax_main.inset_axes([0.65, 0.55, 0.3, 0.4])

for i, eta in enumerate(etas):
    mdot_t = mdot_1myr * (t_myr)**(-eta)
    axins.plot(np.log10(t_myr), np.log10(mdot_t), linewidth=1.5, color=colors[i])

# Calcular el minimo valor alcanzado en t=10 Myr (x=1) para y-lim
min_y = -8.0 - max(etas) * 1.0 
min_y_rounded = np.floor(min_y * 2) / 2

# Limites del Zoom: x entre 0 y 1
axins.set_xlim(0, 1)
axins.set_ylim(min_y_rounded, -8)

# Ticks en pasos de 0.5
ticks_y = np.arange(min_y_rounded, -7.9, 0.5)
axins.set_yticks(ticks_y)

axins.tick_params(axis='both', which='major', labelsize=12)
axins.axvline(0.0, color='gray', linestyle='--', alpha=0.7)
axins.axhline(-8.0, color='gray', linestyle=':', alpha=0.7)
axins.grid(True, which="both", ls="--", alpha=0.5)

# Dibujar el cuadro con líneas que conectan el zoom al gráfico principal
ax_main.indicate_inset_zoom(axins, edgecolor="black")

plt.tight_layout()
plot_path = "hartmann_mdot_evolution.png"
plt.savefig(plot_path, dpi=300)
print(f"Saved plot to {plot_path}")

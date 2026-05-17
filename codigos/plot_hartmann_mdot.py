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
    ax_main.plot(np.log10(t_myr), np.log10(mdot_t), linewidth=2.0, color=colors[i])

# Líneas de referencia
ax_main.axvline(0.0, color='gray', linestyle='--', alpha=0.7)
ax_main.axhline(-8.0, color='gray', linestyle=':', alpha=0.7)

ax_main.set_xlabel(r"$\log_{10}(t \,[\mathrm{Myr}])$", fontsize=12)
ax_main.set_ylabel(r"$\log_{10}(\dot{M}_{\mathrm{gas}} \,[M_\odot/\mathrm{yr}])$", fontsize=12)
ax_main.set_title("Evolución Analítica de la Tasa de Acreción (T Tauri stars)\n" + r"$\dot{M}(t) = 10^{-8} (t / 1\,\mathrm{Myr})^{-\eta}\; M_\odot\,\mathrm{yr}^{-1}$", fontsize=14)
ax_main.grid(True, which="both", ls="--", alpha=0.7)

# Agregar textos de los etas en el lado izquierdo (donde eta=2.8 está arriba)
x_text = -1.5

# Arriba: eta = 2.8
y_top = -8.0 - 2.8 * x_text
ax_main.text(x_text, y_top + 0.3, r'$\eta = 2.8$', color=colors[-1], fontsize=12, fontweight='bold', ha='center', va='bottom')

# Abajo: eta = 1.5
y_bottom = -8.0 - 1.5 * x_text
ax_main.text(x_text, y_bottom - 0.3, r'$\eta = 1.5$', color=colors[0], fontsize=12, fontweight='bold', ha='center', va='top')

# Elimino la leyenda por completo

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

axins.axvline(0.0, color='gray', linestyle='--', alpha=0.7)
axins.axhline(-8.0, color='gray', linestyle=':', alpha=0.7)

axins.grid(True, which="both", ls="--", alpha=0.5)

# Dibujar el cuadro con líneas que conectan el zoom al gráfico principal
ax_main.indicate_inset_zoom(axins, edgecolor="black")

plt.tight_layout()
plot_path = "hartmann_mdot_evolution.png"
plt.savefig(plot_path, dpi=300)
print(f"Saved plot to {plot_path}")

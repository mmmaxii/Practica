import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import os

# ============================================================
# 0. ESTILO DE PLOT (Paper)
# ============================================================
plt.rcParams.update({
    "text.usetex":          False,
    "font.family":          "serif",
    "font.size":            14,
    "axes.labelsize":       16,
    "axes.titlesize":       16,
    "xtick.labelsize":      13,
    "ytick.labelsize":      13,
    "xtick.direction":      "in",
    "ytick.direction":      "in",
    "xtick.top":            True,
    "ytick.right":          True,
    "xtick.minor.visible":  True,
    "ytick.minor.visible":  True,
    "figure.dpi":           150,
})

_tick_kw = dict(
    which="major", direction="in", top=True, right=True,
    length=8, width=1.0, labelsize=13,
)
_tick_kw_minor = dict(
    which="minor", direction="in", top=True, right=True,
    length=4, width=0.7,
)

# ============================================================
# 1. GENERACIÓN DE LA MALLA (GRID) PARA EL GRADIENTE
# ============================================================
# Rango temporal de 100 yr a 10^7 yr
t_yr = np.logspace(2, 7, 500)
t_myr = t_yr / 1e6
log_t = np.log10(t_myr)

# Malla fina de etas para generar el gradiente continuo
etas = np.linspace(1.5, 2.8, 500)
mdot_1myr = 1e-8 # M_sun / yr

# Crear matrices 2D para X (tiempo), Y (mdot) y Z (eta)
T_mat, Eta_mat = np.meshgrid(t_myr, etas)
Log_T_mat = np.log10(T_mat)

# Ecuación de Hartmann (Y es variable dependiente de X y Z)
Mdot_mat = mdot_1myr * (T_mat)**(-Eta_mat)
Log_Mdot_mat = np.log10(Mdot_mat)

# ============================================================
# 2. CREACIÓN DE LA FIGURA
# ============================================================
fig, ax_main = plt.subplots(figsize=(11, 7))

# Colormap para el gradiente (magma o inferno quedan perfectos para densidad)
cmap = plt.colormaps["magma"]

# Plot principal: pcolormesh transforma la grilla X, Y en un gradiente Z suave
mesh_main = ax_main.pcolormesh(Log_T_mat, Log_Mdot_mat, Eta_mat, 
                               cmap=cmap, shading='gouraud', rasterized=True)

# Bordes superior e inferior (las curvas límite)
ax_main.plot(log_t, np.log10(mdot_1myr * t_myr**(-1.5)), color='k', lw=1.8, ls='--', label=r'Límites de $\eta$ (1.5, 2.8)')
ax_main.plot(log_t, np.log10(mdot_1myr * t_myr**(-2.8)), color='k', lw=1.8, ls='--')

# Líneas de referencia
ax_main.axvline(0.0, color='gray', linestyle=':', lw=1.5, alpha=0.8, label='1 Myr')
ax_main.axhline(-8.0, color='gray', linestyle=':', lw=1.5, alpha=0.8)

ax_main.set_xlabel(r"$\log_{10}(t \,[\mathrm{Myr}])$")
ax_main.set_ylabel(r"$\log_{10}(\dot{M}_{\mathrm{gas}} \,[M_\odot/\mathrm{yr}])$")
ax_main.set_title("Evolución de Acreción en Estrellas T Tauri (Gradiente de $\eta$)")
ax_main.tick_params(**_tick_kw)
ax_main.tick_params(**_tick_kw_minor)
ax_main.legend(loc="lower left", fontsize=12, framealpha=0.9, edgecolor="none")

# ============================================================
# 3. ZOOM (Inset Axes)
# ============================================================
axins = ax_main.inset_axes([0.62, 0.52, 0.35, 0.45])

mesh_ins = axins.pcolormesh(Log_T_mat, Log_Mdot_mat, Eta_mat, 
                            cmap=cmap, shading='gouraud', rasterized=True)

# Bordes en el zoom
axins.plot(log_t, np.log10(mdot_1myr * t_myr**(-1.5)), color='k', lw=1.5, ls='--')
axins.plot(log_t, np.log10(mdot_1myr * t_myr**(-2.8)), color='k', lw=1.5, ls='--')

# Calcular el minimo valor alcanzado en t=10 Myr (x=1) para y-lim
min_y = -8.0 - max(etas) * 1.0 
min_y_rounded = np.floor(min_y * 2) / 2

# Limites del Zoom: x entre 0 y 1
axins.set_xlim(0, 1)
axins.set_ylim(min_y_rounded, -8)

# Ticks en pasos de 0.5
ticks_y = np.arange(min_y_rounded, -7.9, 0.5)
axins.set_yticks(ticks_y)

axins.tick_params(**_tick_kw)
axins.tick_params(**_tick_kw_minor)
axins.axvline(0.0, color='gray', linestyle=':', lw=1.5, alpha=0.8)
axins.axhline(-8.0, color='gray', linestyle=':', lw=1.5, alpha=0.8)

# Dibujar el cuadro con líneas que conectan el zoom al gráfico principal
ax_main.indicate_inset_zoom(axins, edgecolor="black", alpha=0.7)

# ============================================================
# 4. COLORBAR
# ============================================================
cbar = fig.colorbar(mesh_main, ax=ax_main, pad=0.02)
cbar.set_label(r"Parámetro $\eta$", fontsize=16, labelpad=10)
cbar.ax.tick_params(labelsize=13)

# Guardar
plot_path = r"c:\astro\Codigos practica + docs + papers\codigos\hartmann_mdot_gradient.png"
pdf_path = plot_path.replace(".png", ".pdf")

plt.savefig(plot_path, dpi=200, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
print(f"\n[OK] Guardado en:\n  {plot_path}\n  {pdf_path}")

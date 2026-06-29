import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import matplotlib.ticker as ticker
import os

# ============================================================
# 0. ESTILO DE PLOT (Paper, carta, 1 columna)
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
    "lines.linewidth":      2.0,
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
# 1. CARGA DE DATOS (Oka et al. 2011)
# ============================================================
csv_path = r"c:\astro\Codigos practica + docs + papers\codigos\Oka et al 2011 fig 7\oka_et_al_2011_data_corrected.csv"

mdot_list = []
rsnow_list = []
with open(csv_path, 'r') as f:
    lines = f.readlines()
    for line in lines[1:]: # Saltar encabezado
        line = line.strip()
        if not line: continue
        parts = line.split(';')
        if len(parts) == 2:
            mdot_list.append(float(parts[0]))
            rsnow_list.append(float(parts[1]))

mdot_raw = np.array(mdot_list)
rsnow_raw = np.array(rsnow_list)

sort_idx = np.argsort(mdot_raw)
mdot_sorted = mdot_raw[sort_idx]
rsnow_sorted = rsnow_raw[sort_idx]

mdot_unique, unique_indices = np.unique(mdot_sorted, return_index=True)
rsnow_unique = rsnow_sorted[unique_indices]

# ============================================================
# 2. INTERPOLACIÓN FÍSICA (Log-Log)
# ============================================================
log_mdot = np.log10(mdot_unique)
log_rsnow = np.log10(rsnow_unique)

interp_log_rsnow = interp1d(
    log_mdot, 
    log_rsnow, 
    kind='linear', 
    bounds_error=False, 
    fill_value="extrapolate"
)

def r_snow_from_mdot(mdot_val):
    log_m = np.log10(mdot_val)
    log_r = interp_log_rsnow(log_m)
    return 10**log_r

# ============================================================
# 3. EVOLUCIÓN TEMPORAL (Hartmann et al.)
# ============================================================
def mdot_time(t_myr, eta=1.5):
    mdot_1myr = 1e-8
    return mdot_1myr * (t_myr)**(-eta)

def r_snow_time(t_myr, eta=1.5):
    m = mdot_time(t_myr, eta)
    return r_snow_from_mdot(m)

# ============================================================
# 4. GRÁFICOS Y DIAGNÓSTICOS
# ============================================================

fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)
ax1, ax2, ax3 = axes

eta_plot = 1.5
color_main = '#d95f02'  # Naranja paper
color_raw = '#1f78b4'   # Azul paper

# --- PANEL 1: Evolución de Mdot(t) ---
t_eval = np.logspace(-1, 1, 500) # 0.1 a 10 Myr
log_t_eval = np.log10(t_eval)

m_t = mdot_time(t_eval, eta=eta_plot)
log_m_t = np.log10(m_t)

ax1.plot(log_t_eval, log_m_t, color=color_main, linewidth=2.5, label=fr'$\eta = {eta_plot}$')
ax1.set_xlim(-1, 1)
# 5 ticks para log(t)
ax1.set_xticks([-1, -0.5, 0, 0.5, 1])
ax1.set_xlabel(r'$\log_{10}(t\ [\mathrm{Myr}])$')
ax1.set_ylabel(r'$\log_{10}(\dot{M}\ [M_\odot\ \mathrm{yr}^{-1}])$')
ax1.legend(loc="upper right", framealpha=0.85, edgecolor="none", fontsize=11)

# --- PANEL 2: Verificación de la Interpolación Oka et al ---
m_test = np.logspace(-12, -7, 500)
r_test = r_snow_from_mdot(m_test)

log_mdot_raw = np.log10(mdot_raw)
log_m_test = np.log10(m_test)

ax2.plot(log_mdot_raw, rsnow_raw, color=color_raw, linestyle='-', linewidth=3, alpha=0.5, label='Oka et al. 2011')
ax2.plot(log_m_test, r_test, color=color_main, linestyle='--', linewidth=2.5, label='Interpolación')

ax2.set_xlim(-7, -12) # Invertido como en el paper

ax2.set_yscale('log')
ax2.set_ylim(0.1, 10)
ax2.set_yticks([0.1, 1, 10])
ax2.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

ax2.set_xlabel(r'$\log_{10}(\dot{M}\ [M_\odot\ \mathrm{yr}^{-1}])$')
ax2.set_ylabel(r'$r_{\mathrm{snow}}\ [\mathrm{AU}]$')

# Linea vertical en log10(Mdot(t=10 Myr))
m_10myr = mdot_time(10, eta=eta_plot)
log_m_10myr = np.log10(m_10myr)
ax2.axvline(log_m_10myr, color='k', linestyle=':', linewidth=1.5, alpha=0.6, label=r'$\dot{M}(10\,\mathrm{Myr})$')

ax2.legend(loc="lower right", framealpha=0.85, edgecolor="none", fontsize=11)

# --- PANEL 3: Evolución de r_snow(t) ---
r_t = r_snow_time(t_eval, eta=eta_plot)

ax3.plot(log_t_eval, r_t, color=color_main, linewidth=2.5, label=fr'$\eta = {eta_plot}$')

ax3.set_xlim(-1, 1)
# 5 ticks para log(t)
ax3.set_xticks([-1, -0.5, 0, 0.5, 1])

ax3.set_yscale('log')
ax3.set_ylim(0.1, 10)
ax3.set_yticks([0.1, 1, 10])
ax3.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

ax3.set_xlabel(r'$\log_{10}(t\ [\mathrm{Myr}])$')
ax3.set_ylabel(r'$r_{\mathrm{snow}}\ [\mathrm{AU}]$')

# Encontrar tiempo donde r_snow = 1 AU
from scipy.optimize import brentq
def root_func(t):
    return r_snow_time(t, eta=eta_plot) - 1.0

try:
    t_1au = brentq(root_func, 0.1, 10.0)
    log_t_1au = np.log10(t_1au)
    
    # Punto negro en el cruce
    ax3.plot(log_t_1au, 1.0, 'ko', markersize=6, zorder=5, label=f'{t_1au:.2f} Myr')
    
    # Lineas punteadas negras para guiar el ojo
    ax3.axvline(log_t_1au, color='k', linestyle=':', linewidth=1.5, alpha=0.6)
    ax3.axhline(1.0, color='k', linestyle=':', linewidth=1.5, alpha=0.6)
except ValueError:
    pass

ax3.legend(loc="lower left", framealpha=0.85, edgecolor="none", fontsize=11)

# Ajustes generales de ticks
for ax in [ax1, ax2, ax3]:
    ax.tick_params(**_tick_kw)
    ax.tick_params(**_tick_kw_minor)
    ax.grid(True, which='both', linestyle=':', alpha=0.3, color="gray")

plot_path = r"c:\astro\Codigos practica + docs + papers\codigos\dynamic_snowline_oka_interpolation.png"
pdf_path = plot_path.replace(".png", ".pdf")

plt.savefig(plot_path, dpi=200, bbox_inches="tight")
plt.savefig(pdf_path, bbox_inches="tight")
print(f"[OK] Guardado en:\n  {plot_path}\n  {pdf_path}")

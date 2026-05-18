import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ============================================================
# DATOS DIGITIZADOS DE OKA ET AL. 2011 (CORREGIDOS)
# ============================================================

csv_path = r"c:\astro\Codigos practica + docs + papers\codigos\Oka et al 2011 fig 7\oka_et_al_2011_data_corrected.csv"

# Leer el CSV manualmente
mdot_list = []
rsnow_list = []
with open(csv_path, 'r') as f:
    lines = f.readlines()
    for line in lines[1:]: # Saltar encabezado
        line = line.strip()
        if not line:
            continue
        parts = line.split(';')
        if len(parts) == 2:
            mdot_list.append(float(parts[0]))
            rsnow_list.append(float(parts[1]))

mdot_all = np.array(mdot_list)
rsnow_all = np.array(rsnow_list)

# Filtrar para hacer el fit solo en el régimen dominado por viscosidad (mdot > 1e-9)
mask = mdot_all > 1e-9
mdot = mdot_all[mask]
rsnow = rsnow_all[mask]

# ============================================================
# MODELO POWER LAW
# r_snow = A * (mdot / 1e-8)^p
# ============================================================

def rsnow_model(mdot, A, p):
    return A * (mdot / 1e-8)**p

# ============================================================
# FIT DIRECTO
# ============================================================

popt, pcov = curve_fit(
    rsnow_model,
    mdot,
    rsnow,
    p0=[1.0, 0.4]
)

A_fit, p_fit = popt

# Errores 1-sigma
A_err, p_err = np.sqrt(np.diag(pcov))

print("\n===== FIT RESULTS =====")
print(f"A = {A_fit:.4f} ± {A_err:.4f}")
print(f"p = {p_fit:.4f} ± {p_err:.4f}")

# ============================================================
# FIT EN ESPACIO LOG-LOG
# ============================================================

log_mdot = np.log10(mdot)
log_rsnow = np.log10(rsnow)

coeffs = np.polyfit(log_mdot, log_rsnow, 1)

p_log = coeffs[0]
intercept = coeffs[1]

# IMPORTANTE:
# porque usamos mdot directamente y no mdot/1e-8
# log(R) = p_log * log(M) + intercept
# Si R = A * (M / 1e-8)^p_log -> log(R) = log(A) + p_log * (log(M) - (-8))
# log(A) + p_log * log(M) + 8 * p_log = intercept + p_log * log(M)
# log(A) + 8 * p_log = intercept => log(A) = intercept - 8 * p_log
# Por lo tanto, A = 10**(intercept - (-8) * p_log) = 10**(intercept + 8 * p_log)
A_log = 10**(intercept + 8 * p_log)

print("\n===== LOG-LOG FIT =====")
print(f"A_log = {A_log:.4f}")
print(f"p_log = {p_log:.4f}")

# ============================================================
# GENERAR CURVA SUAVE
# ============================================================

mdot_fit = np.logspace(
    np.log10(mdot.min()/2),
    np.log10(mdot.max()*2),
    300
)

rs_fit = rsnow_model(mdot_fit, A_fit, p_fit)

# ============================================================
# PLOT
# ============================================================

fig, ax = plt.subplots(figsize=(8,6))

# Dibujar todos los puntos de la data corregida con una linea semi-transparente
ax.plot(
    mdot_all,
    rsnow_all,
    'k-',
    linewidth=4,
    alpha=0.4,
    label="Oka et al. 2011 data"
)

ax.plot(
    mdot_fit,
    rs_fit,
    'r--',
    linewidth=2,
    label=(
        rf"$r_{{snow}} = {A_fit:.2f}"
        rf"\left(\frac{{\dot{{M}}}}{{10^{{-8}}}}\right)^{{{p_fit:.2f}}}$"
    )
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(1e-7, 1e-12)
ax.set_ylim(0.1, 10)

ax.set_xlabel(r"$\dot{M}\ [M_\odot\ yr^{-1}]$", fontsize=14)
ax.set_ylabel(r"$r_{snow}\ [AU]$", fontsize=14)

ax.legend(fontsize=12)

# Estilo de los ticks
ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=12)

plt.tight_layout()
plot_out = r"c:\astro\Codigos practica + docs + papers\codigos\oka_fit_results.png"
plt.savefig(plot_out, dpi=300)
print(f"\nGráfico guardado en: {plot_out}")

# ============================================================
# RESIDUOS
# ============================================================

rs_pred = rsnow_model(mdot, A_fit, p_fit)

residuals = (rsnow - rs_pred) / rsnow

print("\nMean absolute residual:")
print(f"{100*np.mean(np.abs(residuals)):.2f}%")

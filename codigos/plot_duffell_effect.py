"""
plot_duffell_effect.py
======================
Panel 1×3 — Efecto Duffell leído directo desde los HDF5 de simulación.

Datos:
  - 1 caso smooth:  vf_10ms/smooth/run_smooth_a0.001_v10/data0000.hdf5
  - N casos general: vf_10ms/general/run_r10.0_m{M}_a0.001/data0000.hdf5
    (todos los M_gap disponibles, alpha=0.001, r_gap=10 AU, t=0)

Columnas:
  1) alpha(r)
  2) Sigma_gas(r) y Sigma_dust(r)  [polvo = suma sobre bins]
  3) eta(r)

Colormap: inferno  (oscuro = smooth/menor M, amarillo = mayor M)
"""

import os
import re
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import dustpy.constants as c

# ─────────────────────────────────────────────────────────────────────────────
# Rutas base
# ─────────────────────────────────────────────────────────────────────────────

BASE     = "data/runs/vf_10ms"
ALPHA    = "0.001"
R_GAP    = "10.0"
SNAP     = "data0059.hdf5"   # ~ 95,500 yr  (snap mas cercano a 100k yr en log-space)

smooth_path = os.path.join(
    BASE, "smooth", f"run_smooth_a{ALPHA}_v10", SNAP
)

general_dir = os.path.join(BASE, "general")

# ─────────────────────────────────────────────────────────────────────────────
# Descubrir todos los casos general con r=10 AU y a=0.001
# ─────────────────────────────────────────────────────────────────────────────

pattern = re.compile(
    rf"run_r{re.escape(R_GAP)}_m([\d.]+)_a{re.escape(ALPHA)}$"
)

general_cases = []   # list of (M_Mjup, hdf5_path)

for name in os.listdir(general_dir):
    m = pattern.match(name)
    if m:
        M_val  = float(m.group(1))
        hpath  = os.path.join(general_dir, name, SNAP)
        if os.path.isfile(hpath):
            general_cases.append((M_val, hpath))

general_cases.sort(key=lambda x: x[0])   # orden creciente en M_gap

print(f"[INFO] Smooth:   {smooth_path}")
print(f"[INFO] General cases encontrados (r={R_GAP} AU, alpha={ALPHA}):")
for M, p in general_cases:
    print(f"         M = {M:5.2f} M_Jup  ->  {os.path.basename(os.path.dirname(p))}")

# ─────────────────────────────────────────────────────────────────────────────
# Función lectora
# ─────────────────────────────────────────────────────────────────────────────

def read_snapshot(hdf5_path):
    """Lee grid/r, gas/alpha, gas/Sigma, dust/Sigma_total, gas/eta."""
    with h5py.File(hdf5_path, "r") as f:
        r          = f["grid/r"][:]                    # [cm]
        alpha      = f["gas/alpha"][:]
        sigma_gas  = f["gas/Sigma"][:]                 # [g/cm²]
        sigma_dust = f["dust/Sigma"][:].sum(axis=1)   # suma sobre bins de masa
        eta        = f["gas/eta"][:]
        t          = float(f["t"][()])                 # [s]
    r_au = r / c.au
    return dict(r_au=r_au, alpha=alpha,
                sigma_gas=sigma_gas, sigma_dust=sigma_dust, eta=eta, t=t)

# ─────────────────────────────────────────────────────────────────────────────
# Leer todos los perfiles
# ─────────────────────────────────────────────────────────────────────────────

smooth_data = read_snapshot(smooth_path)

cases = []   # list of (label, data_dict)
cases.append((r"smooth  ($M_{\rm gap}=0$)", smooth_data))

for M_val, hpath in general_cases:
    data  = read_snapshot(hpath)
    label = rf"$M_{{\rm gap}} = {M_val}\,M_{{\rm Jup}}$"
    cases.append((label, data))

N = len(cases)
print(f"\n[INFO] Total escenarios: {N}  (1 smooth + {N-1} general)")

# ─────────────────────────────────────────────────────────────────────────────
# Colormap inferno
# ─────────────────────────────────────────────────────────────────────────────

# Colormap inferno truncado: oscuro -> naranja (evita amarillo)
# Muestreamos solo hasta el 82% del rango para quedarnos en naranja
_cmap_base = plt.colormaps["inferno"]
cmap       = mcolors.LinearSegmentedColormap.from_list(
    "inferno_orange",
    _cmap_base(np.linspace(0.0, 0.82, 256))
)
norm   = mcolors.Normalize(vmin=0, vmax=N - 1)
colors = [cmap(norm(i)) for i in range(N)]

# ─────────────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────────────

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
fig, axes = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

ax1, ax2, ax3 = axes

# Ticks estilo paper: major largos y gruesos, minor visibles
_tick_kw = dict(
    which="major", direction="in", top=True, right=True,
    length=8, width=1.0, labelsize=13,
)
_tick_kw_minor = dict(
    which="minor", direction="in", top=True, right=True,
    length=4, width=0.7,
)

gap_vline_kw = dict(color="silver", lw=1.0, ls="--", zorder=0)

# ── Panel 1: alpha ───────────────────────────────────────────────────
for i, (lbl, d) in enumerate(cases):
    lw = 1.8
    ax1.plot(d["r_au"], d["alpha"], color=colors[i], lw=lw, zorder=i + 1)

ax1.axvline(float(R_GAP), **gap_vline_kw)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel(r"$r$ [AU]")
ax1.set_ylabel(r"$\alpha$")
ax1.set_xlim(
    smooth_data["r_au"].min(),
    smooth_data["r_au"].max()
)
ax1.tick_params(**_tick_kw)
ax1.tick_params(**_tick_kw_minor)

# ── Panel 2: Sigma_gas y Sigma_dust ─────────────────────────────────────────
# ── Panel 2: Sigma_gas / Sigma_smooth  y  Sigma_dust / Sigma_smooth ─────────
# Interpolar smooth sobre cada grilla (pueden tener distinto Nr)
r_sm      = smooth_data["r_au"]
sg_sm     = smooth_data["sigma_gas"]
sd_sm     = smooth_data["sigma_dust"]

for i, (lbl, d) in enumerate(cases):
    lw = 2.5 if i == 0 else 1.8
    r_i = d["r_au"]

    # smooth interpolado en la grilla de este caso
    sg_ref = np.interp(r_i, r_sm, sg_sm)
    sd_ref = np.interp(r_i, r_sm, sd_sm)

    ratio_g = d["sigma_gas"]  / sg_ref
    ratio_d = d["sigma_dust"] / sd_ref

    ax2.plot(r_i, ratio_g, color=colors[i], lw=1.8, ls="-")
    ax2.plot(r_i, ratio_d, color=colors[i], lw=1.8, ls=":")

ax2.axhline(1.0, color="k", lw=0.8, ls="-", alpha=0.35, zorder=0)
ax2.axvline(float(R_GAP), **gap_vline_kw)

legend_lines = [
    Line2D([0], [0], color="dimgray", lw=2, ls="-",
           label=r"$\Sigma_{\rm gas}$"),
    Line2D([0], [0], color="dimgray", lw=2, ls=":",
           label=r"$\Sigma_{\rm dust}$"),
]
ax2.legend(handles=legend_lines, fontsize=9, loc="lower left",
           framealpha=0.85, edgecolor="none")


ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel(r"$r$ [AU]")
ax2.set_ylabel(r"$\Sigma\,/\,\Sigma_{\rm smooth}$")
ax2.set_xlim(
    smooth_data["r_au"].min(),
    smooth_data["r_au"].max()
)
ax2.set_ylim(1e-3, 1e2)
ax2.tick_params(**_tick_kw)
ax2.tick_params(**_tick_kw_minor)

# ── Panel 3: eta(r) ──────────────────────────────────────────────────────────
for i, (lbl, d) in enumerate(cases):
    ax3.plot(d["r_au"], d["eta"], color=colors[i], lw=1.8, zorder=i + 1)

ax3.axvline(float(R_GAP), **gap_vline_kw)
ax3.axhline(0.0, color="k", lw=0.7, ls="-", alpha=0.35, zorder=0)

# Sombrear región eta < 0 del caso de mayor M_gap
eta_max = cases[-1][1]["eta"]
r_max   = cases[-1][1]["r_au"]
eta_neg = np.where(eta_max < 0, eta_max, 0.0)
if np.any(eta_max < 0):
    ax3.fill_between(r_max, 0, eta_neg,
                     alpha=0.15, color=colors[-1], label="_nolegend_")
    idx_min = int(np.argmin(eta_max))
    ax3.annotate(
        r"$\eta < 0$" "\n" r"trampa de pebbles",
        xy=(r_max[idx_min], eta_max[idx_min]),
        fontsize=8.5, color=colors[-1], ha="left",
        xytext=(r_max[idx_min] * 1.8, eta_max[idx_min] * 0.6),
        arrowprops=dict(arrowstyle="->", color=colors[-1], lw=0.9),
    )

ax3.set_xscale("log")
ax3.set_xlabel(r"$r$ [AU]")
ax3.set_ylabel(r"$\eta$")
ax3.set_xlim(
    smooth_data["r_au"].min(),
    smooth_data["r_au"].max()
)
ax3.set_ylim(-0.075, 0.075)
ax3.tick_params(**_tick_kw)
ax3.tick_params(**_tick_kw_minor)

# ── Colorbar con punta triangular (estilo paper) ─────────────────────────────
sm = cm.ScalarMappable(
    cmap=cmap,
    norm=mcolors.Normalize(vmin=0, vmax=N - 1)
)
sm.set_array([])
cbar = fig.colorbar(
    sm, ax=axes.tolist(),
    orientation="vertical",
    shrink=0.92, pad=0.015,
    extend="max",          # punta triangular en el extremo superior
)
cbar.set_label(
    r"$M_{\rm gap}\;[M_{\rm Jup}]$",
    fontsize=15, labelpad=10
)
cbar.set_ticks(list(range(N)))
cbar.set_ticklabels(
    ["smooth"] + [str(m) for m, _ in general_cases],
    fontsize=12
)
cbar.ax.tick_params(length=5, width=0.8)

# ── Guardar ──────────────────────────────────────────────────────────────────
outdir = "data/benchmarks"
os.makedirs(outdir, exist_ok=True)

pdf_out = os.path.join(outdir, "duffell_effect_panel.pdf")
png_out = os.path.join(outdir, "duffell_effect_panel.png")

plt.savefig(pdf_out, bbox_inches="tight")
plt.savefig(png_out, dpi=200, bbox_inches="tight")
print(f"\n[OK] Guardado en:\n  {pdf_out}\n  {png_out}")
plt.show()

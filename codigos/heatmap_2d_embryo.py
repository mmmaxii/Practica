"""
heatmap_2d_embryo.py — Mapa 2D: posición del gap × posición del embrión
==================================================================================

Para cada run de tipo single_<masa>_<pos>au en dataroot, corre PA3 con múltiples
posiciones de embrión y construye un mapa 2D por masa de planeta:

    Eje X  →  posición del gap planetario  [AU]
    Eje Y  →  posición del embrión receptor [AU]
    Color  →  fracción final de H₂O [%]  o  masa final [M⊕]

Para cada run de tipo single_<masa>_<pos>au en dataroot, corre PA3 con múltiples
posiciones de embrión y construye un mapa 2D por masa de planeta:

    Eje X  →  posición del gap planetario  [AU]
    Eje Y  →  posición del embrión receptor [AU]
    Color  →  fracción final de H₂O  [%]

Genera un subplot por tipo de masa (Super-Earth, Neptune, Saturn, Jupiter, Super-Jupiter)
más un panel "promedio" sobre todas las masas.

Uso
---
    py heatmap_2d_embryo.py
    py heatmap_2d_embryo.py --dataroot data/5myr --savedir figures/5myr --cmap inferno
    py heatmap_2d_embryo.py --quantity mass --cmap viridis

Argumentos
----------
  --dataroot   Carpeta raíz con los runs  (default: data/5myr)
  --savedir    Carpeta de salida          (default: figures/5myr)
  --cmap       Colormap matplotlib        (default: inferno)
  --quantity   'fh2o' (fracción H₂O %) o 'mass' (masa final M⊕) (default: fh2o)
  --embryos    Posiciones del embrión [AU] (default: 1.0 2.0 3.0 5.0 7.0 10.0)
  --t_min_yr   Primer snapshot a usar [yr] (default: 100.0)
  --vmax       Máximo del colorbar         (default: auto)
  --ww         Umbral waterworld [%]       (default: 10.0, solo para fh2o)
"""

import os
import sys
import re
import glob
import argparse

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable

import dustpy.constants as c

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3

# ── Constantes ────────────────────────────────────────────────────────────────
AU      = c.au
YR      = c.year
M_EARTH = c.M_earth

# ── Estética oscura ───────────────────────────────────────────────────────────
BG     = "#0d1117"
FG     = "#e6edf3"
GRID_C = "#21262d"
ACCENT = "#ffd700"

plt.rcParams.update({
    "figure.facecolor":  BG,
    "axes.facecolor":    "#161b22",
    "axes.edgecolor":    GRID_C,
    "axes.labelcolor":   FG,
    "axes.titlecolor":   FG,
    "xtick.color":       FG,
    "ytick.color":       FG,
    "text.color":        FG,
    "savefig.facecolor": BG,
    "savefig.dpi":       180,
    "font.family":       "sans-serif",
    "font.size":         10,
})

# ── Metadatos de masa ─────────────────────────────────────────────────────────
MASS_ORDER  = ["sup_earth", "nep", "sat", "jup", "sup_jup"]
MASS_LABELS = {
    "sup_earth": "Super-Earth\n(5 M⊕)",
    "nep":       "Neptune\n(17 M⊕)",
    "sat":       "Saturn\n(95 M⊕)",
    "jup":       "Jupiter\n(318 M⊕)",
    "sup_jup":   "Super-Jupiter\n(954 M⊕)",
}


# ══════════════════════════════════════════════════════════════════════════════
# Utilidades
# ══════════════════════════════════════════════════════════════════════════════

def parse_single(name):
    """'single_jup_3.0au' → {'mass': 'jup', 'pos': 3.0} o None"""
    m = re.match(r"single_([a-z_]+)_([\d.]+)au", name)
    if not m:
        return None
    raw = m.group(1)
    if "sup_earth" in raw:
        mass = "sup_earth"
    elif "sup_jup" in raw:
        mass = "sup_jup"
    else:
        mass = raw
    return {"mass": mass, "pos": float(m.group(2))}


def collect_single_runs(dataroot):
    """Devuelve lista de dicts con name, path, mass, gap_pos."""
    runs = []
    for entry in sorted(os.listdir(dataroot)):
        path = os.path.join(dataroot, entry)
        if not os.path.isdir(path):
            continue
        if not glob.glob(os.path.join(path, "data*.hdf5")):
            continue
        meta = parse_single(entry)
        if meta is None:
            continue
        runs.append({"name": entry, "path": path,
                     "mass": meta["mass"], "gap_pos": meta["pos"]})
    return runs


def run_pa3_multi(path, embryo_positions_au, t_min_yr=100.0):
    """Corre PA3 para múltiples posiciones de embrión. Devuelve dict {r_au: hist}."""
    pa3  = PebbleAccretionModule3.from_datadir(path, t_min_yr=t_min_yr)
    results = pa3.run_growth(embryo_positions_au)
    return results


def fh2o_final(hist):
    """Fracción de H₂O final [%]."""
    if hist is None or len(hist) == 0:
        return np.nan
    M = hist[-1, 1]
    H = hist[-1, 2]
    return 100.0 * H / (M + 1e-30) if M > 0 else np.nan


def mass_final(hist):
    """Masa final del embrión [M⊕]."""
    if hist is None or len(hist) == 0:
        return np.nan
    return hist[-1, 1] / M_EARTH


def extract_metric(hist, quantity):
    """Extrae la métrica según quantity: 'fh2o' o 'mass'."""
    if quantity == "mass":
        return mass_final(hist)
    return fh2o_final(hist)


# ══════════════════════════════════════════════════════════════════════════════
# Construcción de matrices de datos
# ══════════════════════════════════════════════════════════════════════════════

def build_matrices(runs, embryo_positions, t_min_yr, gap_positions, quantity="fh2o"):
    """
    Devuelve dict: mass_key → Z (n_embryo × n_gap) con la métrica elegida.
    También devuelve 'mean' con el promedio sobre masas.
    """
    matrices = {m: np.full((len(embryo_positions), len(gap_positions)), np.nan)
                for m in MASS_ORDER}

    n_runs = len(runs)
    for i, run in enumerate(runs):
        mass    = run["mass"]
        gap_pos = run["gap_pos"]

        if gap_pos not in gap_positions:
            continue
        gi = gap_positions.index(gap_pos)

        print(f"  [{i+1}/{n_runs}] PA3 → {run['name']} ...")
        try:
            results = run_pa3_multi(run["path"], embryo_positions, t_min_yr)
        except Exception as e:
            print(f"    [!] Error: {e}")
            continue

        for ei, r_emb in enumerate(embryo_positions):
            hist = results.get(r_emb)
            matrices[mass][ei, gi] = extract_metric(hist, quantity)

    # Promedio sobre masas (ignorando NaN)
    stack = np.stack([matrices[m] for m in MASS_ORDER], axis=0)
    matrices["mean"] = np.nanmean(stack, axis=0)

    return matrices


# ══════════════════════════════════════════════════════════════════════════════
# Plot: panel individual
# ══════════════════════════════════════════════════════════════════════════════

def _draw_panel(ax, Z, gap_positions, embryo_positions, cmap_name, vmin, vmax,
                ww_threshold, annotate=True, quantity="fh2o"):
    """Dibuja un panel del heatmap. Devuelve el mappable para la colorbar."""

    # Para 'mass': escala logarítmica opcional
    norm = None
    if quantity == "mass":
        pos_vals = Z[Z > 0]
        if len(pos_vals) > 0:
            norm = mcolors.LogNorm(vmin=max(vmin, 1e-3), vmax=vmax)
        else:
            norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    cmap = matplotlib.colormaps[cmap_name].copy()
    cmap.set_bad("#0d1117")

    im = ax.imshow(
        Z,
        cmap=cmap,
        norm=norm,
        aspect="auto",
        origin="lower",
        extent=[
            -0.5, len(gap_positions) - 0.5,
            -0.5, len(embryo_positions) - 0.5,
        ],
    )

    # Anotaciones en cada celda
    if annotate:
        for ei in range(len(embryo_positions)):
            for gi in range(len(gap_positions)):
                val = Z[ei, gi]
                if np.isnan(val):
                    continue
                # Elegir color de texto según brillo relativo
                if quantity == "mass":
                    # escala log: comparar con rango
                    log_range = np.log10(max(vmax, 1e-3)) - np.log10(max(vmin if vmin > 0 else 1e-3, 1e-6))
                    log_val   = np.log10(max(val, 1e-6)) - np.log10(max(vmin if vmin > 0 else 1e-3, 1e-6))
                    norm_val  = log_val / max(log_range, 1e-6)
                    label_str = f"{val:.2f} M⊕"
                else:
                    norm_val  = (val - 0) / max(vmax, 1e-6)
                    label_str = f"{val:.1f}%"

                txt_color = "#0d1117" if norm_val > 0.55 else FG
                ax.text(gi, ei, label_str,
                        ha="center", va="center",
                        fontsize=8, fontweight="bold", color=txt_color)

                # Borde dorado solo para fh2o waterworlds
                if quantity == "fh2o" and val >= ww_threshold:
                    ax.add_patch(Rectangle(
                        (gi - 0.5, ei - 0.5), 1, 1,
                        fill=False, edgecolor=ACCENT, lw=2.0, zorder=5
                    ))

    # Ticks
    ax.set_xticks(range(len(gap_positions)))
    ax.set_xticklabels([f"{p:.0f}" for p in gap_positions], fontsize=9)
    ax.set_yticks(range(len(embryo_positions)))
    ax.set_yticklabels(
        [f"{r:.1f}" if r != int(r) else f"{int(r)}" for r in embryo_positions],
        fontsize=9
    )

    # Grid entre celdas
    for x in np.arange(-0.5, len(gap_positions), 1):
        ax.axvline(x, color=GRID_C, lw=0.6, zorder=3)
    for y in np.arange(-0.5, len(embryo_positions), 1):
        ax.axhline(y, color=GRID_C, lw=0.6, zorder=3)

    return im


# ══════════════════════════════════════════════════════════════════════════════
# Figura principal
# ══════════════════════════════════════════════════════════════════════════════

def plot_2d_heatmap(matrices, gap_positions, embryo_positions,
                    cmap_name, vmax_user, ww_threshold, savedir, tag="", quantity="fh2o"):
    """
    6 paneles: 5 masas + 1 promedio.
    """
    panels = MASS_ORDER + ["mean"]
    panel_titles = {**MASS_LABELS,
                    "mean": "Promedio\n(todas las masas)"}

    # Escala
    all_vals = np.concatenate([matrices[k].ravel() for k in panels])
    valid    = all_vals[~np.isnan(all_vals) & (all_vals > 0)]
    if quantity == "mass":
        vmin = float(np.nanmin(valid)) if len(valid) > 0 else 1e-3
        vmax = vmax_user if vmax_user > 0 else float(np.nanpercentile(valid, 98)) if len(valid) > 0 else 10.0
        cbar_label = r"$M_{\rm final}$ [$M_\oplus$]"
        fname_base = "heatmap_2d_mass"
    else:
        vmin = 0.0
        vmax = vmax_user if vmax_user > 0 else float(np.nanpercentile(all_vals[~np.isnan(all_vals)], 98))
        vmax = max(vmax, 5.0)
        cbar_label = r"$f_{\rm H_2O}$ final [%]"
        fname_base = "heatmap_2d_embryo"

    n_cols = 3
    n_rows = 2
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(6 * n_cols, 5.5 * n_rows),
                             facecolor=BG)
    axes = axes.ravel()

    shared_im = None
    for idx, key in enumerate(panels):
        ax = axes[idx]
        Z  = matrices[key]
        im = _draw_panel(ax, Z, gap_positions, embryo_positions,
                         cmap_name, vmin, vmax, ww_threshold,
                         annotate=True, quantity=quantity)
        shared_im = im

        title = panel_titles.get(key, key)
        style = {"fontsize": 11, "fontweight": "bold", "color": FG, "pad": 8}
        if key == "mean":
            style["color"] = ACCENT
        ax.set_title(title, **style)
        ax.set_xlabel("Posición del gap [AU]", fontsize=9.5, color=FG)
        ax.set_ylabel("Posición del embrión [AU]", fontsize=9.5, color=FG)
        ax.tick_params(colors=FG)

    # Ocultar eje extra si n_panels < n_cols*n_rows
    for idx in range(len(panels), len(axes)):
        axes[idx].set_visible(False)

    # Colorbar compartida
    fig.subplots_adjust(right=0.88, hspace=0.45, wspace=0.35)
    cbar_ax = fig.add_axes([0.91, 0.12, 0.018, 0.76])
    cb = fig.colorbar(shared_im, cax=cbar_ax)
    cb.set_label(cbar_label, color=FG, fontsize=12)
    cb.ax.tick_params(colors=FG, labelsize=10)
    cb.outline.set_edgecolor(GRID_C)

    # Línea waterworld en la colorbar (solo f_H2O)
    if quantity == "fh2o":
        cb.ax.axhline(ww_threshold, color=ACCENT, lw=2, ls="--")
        cb.ax.text(2.5, ww_threshold, f"WW ≥{ww_threshold:.0f}%",
                   color=ACCENT, fontsize=8, va="center")

    # Título global
    dataname = os.path.basename(savedir) if savedir else ""
    fig.suptitle(
        f"Fracción de H₂O final  |  gap pos × embryo pos\n"
        f"Colormap: {cmap_name}  —  datos: {dataname}  —  borde dorado ≥{ww_threshold:.0f}%",
        fontsize=13, fontweight="bold", color=FG, y=1.01
    )

    fname = f"{fname_base}{'_' + tag if tag else ''}.png"
    out   = os.path.join(savedir, fname)
    os.makedirs(savedir, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=180)
    print(f"\n[OK] Heatmap 2D guardado: {out}")
    plt.close(fig)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# Figura compacta: solo el panel "mean" en grande
# ══════════════════════════════════════════════════════════════════════════════

def plot_mean_only(matrices, gap_positions, embryo_positions,
                   cmap_name, vmax_user, ww_threshold, savedir, tag="", quantity="fh2o"):
    """Genera solo el panel de promedio, más grande y limpio."""
    Z    = matrices["mean"]
    valid = Z[~np.isnan(Z) & (Z > 0)]
    if quantity == "mass":
        vmin = float(np.nanmin(valid)) if len(valid) > 0 else 1e-3
        vmax = vmax_user if vmax_user > 0 else float(np.nanpercentile(valid, 98)) if len(valid) > 0 else 10.0
        cbar_label = r"$M_{\rm final}$ [$M_\oplus$]"
        title_qty  = "Masa final promedio [M\u2295]"
        fname_base = "heatmap_2d_mean_mass"
    else:
        vmin = 0.0
        vmax = vmax_user if vmax_user > 0 else float(np.nanpercentile(Z[~np.isnan(Z)], 98))
        vmax = max(vmax, 5.0)
        cbar_label = r"$f_{\rm H_2O}$ [%]"
        title_qty  = "Fracción de H₂O promedio (%)"
        fname_base = "heatmap_2d_mean"

    fig, ax = plt.subplots(figsize=(9, 6), facecolor=BG)
    im = _draw_panel(ax, Z, gap_positions, embryo_positions,
                     cmap_name, vmin, vmax, ww_threshold,
                     annotate=True, quantity=quantity)

    ax.set_title(title_qty + "\n(promediado sobre todas las masas de planeta)",
                 fontsize=13, fontweight="bold", color=FG, pad=10)
    ax.set_xlabel("Posición del gap planetario [AU]", fontsize=11, color=FG)
    ax.set_ylabel("Posición del embrión [AU]", fontsize=11, color=FG)
    ax.tick_params(colors=FG)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="4%", pad=0.12)
    cb  = fig.colorbar(im, cax=cax)
    cb.set_label(cbar_label, color=FG, fontsize=11)
    cb.ax.tick_params(colors=FG)
    cb.outline.set_edgecolor(GRID_C)
    if quantity == "fh2o":
        cb.ax.axhline(ww_threshold, color=ACCENT, lw=1.8, ls="--")

    fig.tight_layout()
    fname = f"{fname_base}{'_' + tag if tag else ''}.png"
    out   = os.path.join(savedir, fname)
    os.makedirs(savedir, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=180)
    print(f"[OK] Panel promedio guardado: {out}")
    plt.close(fig)
    return out


# ══════════════════════════════════════════════════════════════════════════════
# Figura comparativa multi-dataset (ej. 1 Myr vs 5 Myr)
# ══════════════════════════════════════════════════════════════════════════════

def plot_comparison_datasets(datasets, embryo_positions, t_min_yr,
                             cmap_name, vmax_user, ww_threshold,
                             savedir, quantity="fh2o"):
    """
    datasets: list de dicts {'dataroot': str, 'label': str}
    Genera una figura con N columnas (una por dataset) mostrando el panel 'mean'.
    También genera una figura de diferencia (dataset[-1] - dataset[0]) si son 2.
    """
    n_ds = len(datasets)
    all_matrices = []
    common_gaps  = None

    for ds in datasets:
        print(f"\n  >> Dataset: {ds['label']}  ({ds['dataroot']})")
        runs = collect_single_runs(ds["dataroot"])
        if not runs:
            print(f"  [!] Sin runs en {ds['dataroot']} — saltando")
            all_matrices.append(None)
            continue
        gap_positions = sorted({r["gap_pos"] for r in runs})
        if common_gaps is None:
            common_gaps = gap_positions
        mat = build_matrices(runs, embryo_positions, t_min_yr, gap_positions, quantity=quantity)
        mat["_gaps"] = gap_positions
        all_matrices.append(mat)

    if common_gaps is None:
        print("[ERROR] Ningún dataset válido.")
        return

    # ── Escala global compartida ────────────────────────────────────────────
    all_vals = np.concatenate([
        m["mean"].ravel() for m in all_matrices if m is not None
    ])
    valid = all_vals[~np.isnan(all_vals) & (all_vals > 0)]
    if quantity == "mass":
        vmin_g = float(np.nanmin(valid)) if len(valid) > 0 else 1e-3
        vmax_g = vmax_user if vmax_user > 0 else float(np.nanpercentile(valid, 98)) if len(valid) > 0 else 10.0
        cbar_label = r"$M_{\rm final}$ [$M_\oplus$]"
        fname_base = "compare_datasets_mass"
    else:
        vmin_g = 0.0
        vmax_g = vmax_user if vmax_user > 0 else float(np.nanpercentile(all_vals[~np.isnan(all_vals)], 98))
        vmax_g = max(vmax_g, 5.0)
        cbar_label = r"$f_{\rm H_2O}$ [%]"
        fname_base = "compare_datasets_fh2o"

    # ── Figura principal: N columnas (panel mean de cada dataset) ───────────
    n_cols_fig = n_ds + (1 if n_ds == 2 else 0)  # +1 para diferencia
    fig, axes = plt.subplots(1, n_cols_fig,
                             figsize=(8 * n_cols_fig, 6.5),
                             facecolor=BG)
    if n_cols_fig == 1:
        axes = [axes]

    shared_im = None
    for col, (ds, mat) in enumerate(zip(datasets, all_matrices)):
        ax = axes[col]
        if mat is None:
            ax.set_visible(False)
            continue
        gaps = mat["_gaps"]
        Z    = mat["mean"]
        im   = _draw_panel(ax, Z, gaps, embryo_positions,
                           cmap_name, vmin_g, vmax_g, ww_threshold,
                           annotate=True, quantity=quantity)
        shared_im = im
        ax.set_title(ds["label"], fontsize=13, fontweight="bold",
                     color=ACCENT, pad=10)
        ax.set_xlabel("Posición del gap [AU]", fontsize=10, color=FG)
        ax.set_ylabel("Posición del embrión [AU]", fontsize=10, color=FG)
        ax.tick_params(colors=FG)

    # ── Panel de diferencia (solo si hay exactamente 2 datasets) ───────────
    if n_ds == 2 and all_matrices[0] is not None and all_matrices[1] is not None:
        ax_diff = axes[2]

        # Alinear a los gaps comunes entre los dos datasets
        gaps0 = all_matrices[0]["_gaps"]
        gaps1 = all_matrices[1]["_gaps"]
        common_gaps_diff = sorted(set(gaps0) & set(gaps1))

        if not common_gaps_diff:
            ax_diff.text(0.5, 0.5, "Sin gaps comunes\nentre datasets",
                         ha="center", va="center", color=FG,
                         fontsize=12, transform=ax_diff.transAxes)
            ax_diff.set_title("Δ  (sin datos)", fontsize=13, color="#78c1ff", pad=10)
        else:
            # Extraer columnas comunes de cada matriz "mean"
            idx0 = [gaps0.index(g) for g in common_gaps_diff]
            idx1 = [gaps1.index(g) for g in common_gaps_diff]
            Z0_common = all_matrices[0]["mean"][:, idx0]
            Z1_common = all_matrices[1]["mean"][:, idx1]
            Z_diff = Z1_common - Z0_common

            vd = np.nanmax(np.abs(Z_diff[~np.isnan(Z_diff)])) * 1.05 if np.any(~np.isnan(Z_diff)) else 1.0
            cmap_div = matplotlib.colormaps["RdBu_r"].copy()
            cmap_div.set_bad("#0d1117")
            norm_div = mcolors.TwoSlopeNorm(vmin=-vd, vcenter=0, vmax=vd)

            im_diff = ax_diff.imshow(
                Z_diff, cmap=cmap_div, norm=norm_div,
                aspect="auto", origin="lower",
                extent=[-0.5, len(common_gaps_diff) - 0.5,
                        -0.5, len(embryo_positions) - 0.5]
            )
            # Anotaciones Δ
            for ei in range(Z_diff.shape[0]):
                for gi in range(Z_diff.shape[1]):
                    val = Z_diff[ei, gi]
                    if np.isnan(val):
                        continue
                    sign = "+" if val >= 0 else ""
                    lbl = f"{sign}{val:.2f}" if quantity == "mass" else f"{sign}{val:.1f}%"
                    norm_v = (val + vd) / (2 * vd + 1e-9)
                    tc = "#0d1117" if 0.35 < norm_v < 0.65 else FG
                    ax_diff.text(gi, ei, lbl, ha="center", va="center",
                                 fontsize=8, fontweight="bold", color=tc)
            ax_diff.set_xticks(range(len(common_gaps_diff)))
            ax_diff.set_xticklabels([f"{g:.0f}" for g in common_gaps_diff], fontsize=9)
            ax_diff.set_yticks(range(len(embryo_positions)))
            ax_diff.set_yticklabels(
                [f"{r:.1f}" if r != int(r) else f"{int(r)}" for r in embryo_positions],
                fontsize=9
            )
            for x in np.arange(-0.5, len(common_gaps_diff), 1):
                ax_diff.axvline(x, color=GRID_C, lw=0.6)
            for y in np.arange(-0.5, len(embryo_positions), 1):
                ax_diff.axhline(y, color=GRID_C, lw=0.6)

            divider = make_axes_locatable(ax_diff)
            cax_diff = divider.append_axes("right", size="4%", pad=0.12)
            cb_diff  = fig.colorbar(im_diff, cax=cax_diff)
            cb_diff.set_label(f"Δ {cbar_label}", color=FG, fontsize=10)
            cb_diff.ax.tick_params(colors=FG)
            cb_diff.outline.set_edgecolor(GRID_C)

        label_b = datasets[1]["label"]
        label_a = datasets[0]["label"]
        ax_diff.set_title(f"Δ  ({label_b} − {label_a})",
                          fontsize=13, fontweight="bold", color="#78c1ff", pad=10)
        ax_diff.set_xlabel("Posición del gap [AU]", fontsize=10, color=FG)
        ax_diff.set_ylabel("Posición del embrión [AU]", fontsize=10, color=FG)
        ax_diff.tick_params(colors=FG)

    # ── Colorbar global (primera/única) ────────────────────────────────────
    if shared_im is not None:
        divider0 = make_axes_locatable(axes[0])
        cax0 = divider0.append_axes("left", size="4%", pad=0.5)
        cb0  = fig.colorbar(shared_im, cax=cax0)
        cb0.set_label(cbar_label, color=FG, fontsize=10)
        cb0.ax.tick_params(colors=FG, labelsize=9)
        cb0.ax.yaxis.set_ticks_position("left")
        cb0.ax.yaxis.set_label_position("left")
        cb0.outline.set_edgecolor(GRID_C)
        if quantity == "fh2o":
            cb0.ax.axhline(ww_threshold, color=ACCENT, lw=1.5, ls="--")

    ds_labels = "  vs  ".join(ds["label"] for ds in datasets)
    fig.suptitle(
        f"Comparación de datasets  —  {ds_labels}\n"
        f"({'Fracción H₂O' if quantity == 'fh2o' else 'Masa final'})  —  promedio sobre masas de gap",
        fontsize=13, fontweight="bold", color=FG, y=1.02
    )
    fig.tight_layout()

    out = os.path.join(savedir, f"{fname_base}.png")
    os.makedirs(savedir, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", dpi=180)
    print(f"\n[OK] Comparación guardada: {out}")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    p = argparse.ArgumentParser(
        description="Heatmap 2D: posición gap × posición embrión → f_H₂O"
    )
    p.add_argument("--dataroot",  default="data/5myr",
                   help="Carpeta raíz con los runs (default: data/5myr)")
    p.add_argument("--datasets",  nargs="+", metavar="RUTA:ETIQUETA",
                   help="Comparar múltiples datasets. Formato: data/1myr:1Myr data/5myr:5Myr")
    p.add_argument("--savedir",   default="figures/comparativa",
                   help="Carpeta de salida (default: figures/comparativa)")
    p.add_argument("--cmap",      default="inferno",
                   help="Colormap matplotlib (default: inferno). Prueba: 'magma', 'viridis', 'plasma'")
    p.add_argument("--quantity",  default="fh2o", choices=["fh2o", "mass"],
                   help="Métrica a visualizar: 'fh2o' (fracción H₂O %%) o 'mass' (masa M⊕) (default: fh2o)")
    p.add_argument("--embryos",   nargs="+", type=float,
                   default=[2.0, 2.5, 3.0, 5.0, 7.0, 10.0],
                   help="Posiciones del embrión [AU] (default: 2.0 2.5 3.0 5.0 7.0 10.0)")
    p.add_argument("--t_min_yr",  type=float, default=100.0,
                   help="Primer snapshot a usar [yr] (default: 100.0)")
    p.add_argument("--vmax",      type=float, default=0.0,
                   help="Máximo del colorbar [%%]. 0 = auto (default: auto)")
    p.add_argument("--ww",        type=float, default=10.0,
                   help="Umbral waterworld [%%] (default: 10.0)")
    p.add_argument("--tag",       default="",
                   help="Sufijo para el nombre del archivo de salida")
    args = p.parse_args()

    print(f"\n{'='*65}")
    print(f"  Heatmap 2D: gap position × embryo position")
    print(f"  quantity : {args.quantity}")
    print(f"  cmap     : {args.cmap}")
    print(f"  savedir  : {args.savedir}")
    print(f"{'='*65}\n")

    # ── Modo comparación multi-dataset ────────────────────────────────────────
    if args.datasets:
        parsed = []
        for item in args.datasets:
            if ":" in item:
                path, label = item.rsplit(":", 1)
            else:
                path  = item
                label = os.path.basename(item.rstrip("/\\"))
            parsed.append({"dataroot": path, "label": label})

        print(f"  Modo comparación: {[d['label'] for d in parsed]}")
        plot_comparison_datasets(
            parsed, args.embryos, args.t_min_yr,
            args.cmap, args.vmax, args.ww,
            args.savedir, quantity=args.quantity
        )
        print(f"\n  Figuras guardadas en {args.savedir}/")
        return

    # ── Modo dataset único ────────────────────────────────────────────────────
    print(f"  dataroot : {args.dataroot}")
    print(f"  embryos  : {args.embryos} AU\n")

    # Recopilar runs
    runs = collect_single_runs(args.dataroot)
    if not runs:
        print(f"[ERROR] No se encontraron runs single_* en {args.dataroot}")
        sys.exit(1)

    gap_positions = sorted({r["gap_pos"] for r in runs})
    print(f"  {len(runs)} runs encontradas  |  gaps en {gap_positions} AU\n")

    # Calcular matrices
    matrices = build_matrices(
        runs, args.embryos, args.t_min_yr, gap_positions, quantity=args.quantity
    )

    # Generar figuras
    plot_2d_heatmap(
        matrices, gap_positions, args.embryos,
        args.cmap, args.vmax, args.ww, args.savedir, args.tag,
        quantity=args.quantity
    )
    plot_mean_only(
        matrices, gap_positions, args.embryos,
        args.cmap, args.vmax, args.ww, args.savedir, args.tag,
        quantity=args.quantity
    )

    print(f"\n  Figuras guardadas en {args.savedir}/")


if __name__ == "__main__":
    main()

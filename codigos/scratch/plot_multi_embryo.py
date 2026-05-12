"""
plot_multi_embryo.py
====================
Pipeline de plots comparativos multi-embrión para todas las runs de t_5e6.

Para cada run calcula PA3 de forma INDEPENDIENTE por posición (sin depleción
cruzada de flujo entre embriones). Se generan 6 figuras:

  Fig 1 — M_final   vs  a [AU]        (perfil radial final)
  Fig 2 — f_H2O_f   vs  a [AU]        (perfil radial H2O final)
  Fig 3 — f_H2O(t)  grilla multi-panel (un panel por posición)
  Fig 4 — M_core(t) grilla multi-panel (un panel por posición)
  Fig 5 — Stacked-bar arquitectura final agrupado por posición
  Fig 6 — Heatmap 2D  runs × posiciones  (color = f_H2O_final)

Codificación visual:
  Color  → run / gap-config  (paleta de 10 colores)
  Marker → posición orbital  (forma geométrica, una por posición)

Uso:
    python plot_multi_embryo.py [--dataroot data_gaps_pipeline/t_5e6]
                                [--savedir  figs_pa3/multi_embryo]
                                [--t_min_yr 100]
                                [--no-show]
"""

from __future__ import annotations

import os
import sys
import glob
import argparse
import warnings
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker   as mticker
import matplotlib.patches  as mpatches
import matplotlib.gridspec as gridspec
import matplotlib.colors   as mcolors
from matplotlib.lines   import Line2D
from matplotlib.patches import Patch

matplotlib.rcParams.update({
    "font.family":      "DejaVu Sans",
    "axes.titlesize":   12,
    "axes.labelsize":   10,
    "xtick.labelsize":  8,
    "ytick.labelsize":  8,
    "axes.grid":        True,
    "grid.alpha":       0.18,
    "grid.linestyle":   "--",
    "figure.dpi":       140,
    "axes.spines.top":  False,
    "axes.spines.right":False,
    "legend.framealpha": 0.25,
    "legend.fontsize":  8,
})

# ── path PA3Py ────────────────────────────────────────────────────────────────
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3

import dustpy.constants as c
M_EARTH = c.M_earth
AU      = c.au
YR      = c.year


# ════════════════════════════════════════════════════════════════════════════
# Diseño visual
# ════════════════════════════════════════════════════════════════════════════

# Colores: uno por run (hasta 10)
_PALETTE = [
    "#4FC3F7",  # azul
    "#FF8A65",  # naranja
    "#81C784",  # verde
    "#CE93D8",  # lila
    "#FFCC80",  # amarillo cálido
    "#80CBC4",  # teal
    "#F48FB1",  # rosa
    "#AED581",  # lima
    "#90CAF9",  # azul pálido
    "#FFAB91",  # melocotón
]

# Marcadores: uno por posición orbital
_MARKERS: dict[float, str] = {
    1.0:  "o",   # círculo
    2.0:  "s",   # cuadrado
    2.5:  "^",   # triángulo arriba
    3.0:  "D",   # diamante
    5.0:  "v",   # triángulo abajo
    10.0: "P",   # plus relleno (filled plus)
}

# Tamaño del marker en scatter / line plots
_MS = 9

# Posiciones orbitales deseadas
_POSITIONS_AU: list[float] = [1.0, 2.0, 2.5, 3.0, 5.0, 10.0]

# Colores de composición (stacked bar, heatmap)
_COMP_COLS = {
    "silicates": "#EF9A9A",
    "CO2":       "#80CBC4",
    "H2O":       "#4FC3F7",
}

# Umbral waterworld
_WW_THRESHOLD = 10.0   # %


# ════════════════════════════════════════════════════════════════════════════
# Utilidades
# ════════════════════════════════════════════════════════════════════════════

def _pretty(name: str) -> str:
    return name.replace("_", " ")


def _collect_runs(dataroot: str) -> list[tuple[str, str]]:
    """Devuelve lista de (nombre, path) con al menos un data*.hdf5."""
    runs = []
    for entry in sorted(os.listdir(dataroot)):
        path = os.path.join(dataroot, entry)
        if not os.path.isdir(path):
            continue
        if glob.glob(os.path.join(path, "data*.hdf5")):
            runs.append((entry, path))
        else:
            print(f"  [skip] {entry}  (sin archivos data*.hdf5)")
    return runs


def _run_independent(run_path: str, positions: list[float], t_min_yr: float) -> dict:
    """
    Para CADA posición orbital corre PA3 de forma INDEPENDIENTE
    (sin depleción cruzada entre embriones).

    Retorna:
        dict  { r_au: np.ndarray shape (Nt, 7) }
        donde cols = [t_s, M_core_g, M_H2O_g, M_CO2_g, M_sil_g, r_snow_AU, M_iso_g]
    """
    pa3 = PebbleAccretionModule3.from_datadir(run_path, t_min_yr=t_min_yr)
    results = {}
    for r_au in positions:
        res = pa3.run_growth([r_au])   # solo este embrión → no hay depleción cruzada
        results[r_au] = res[r_au]
    return results


# ════════════════════════════════════════════════════════════════════════════
# Extracción de métricas finales
# ════════════════════════════════════════════════════════════════════════════

def _final_metrics(hist: np.ndarray) -> dict:
    """Métricas del último snapshot de un historial."""
    if hist is None or len(hist) == 0:
        return {"M_core": np.nan, "f_H2O": np.nan, "f_CO2": np.nan, "f_sil": np.nan}
    row   = hist[-1]
    M     = row[1]
    M_h2o = row[2]
    M_co2 = row[3]
    M_sil = max(0.0, M - M_h2o - M_co2)
    f_h2o = 100.0 * M_h2o / (M + 1e-30)
    f_co2 = 100.0 * M_co2 / (M + 1e-30)
    f_sil = 100.0 * M_sil / (M + 1e-30)
    return {
        "M_core": M / M_EARTH,
        "f_H2O":  f_h2o,
        "f_CO2":  f_co2,
        "f_sil":  f_sil,
    }


def _save(fig, savedir: str, name: str):
    os.makedirs(savedir, exist_ok=True)
    path = os.path.join(savedir, name)
    fig.savefig(path, bbox_inches="tight")
    print(f"  → {path}")


# ════════════════════════════════════════════════════════════════════════════
# Fig 1: M_final  vs  a [AU]
# ════════════════════════════════════════════════════════════════════════════

def fig_mass_radial(all_data: list[dict], savedir: str):
    """
    Perfil radial de masa final.  Un punto por (run, posición).
    Color = run,  Marker = posición.
    """
    fig, ax = plt.subplots(figsize=(10, 5.5))

    for rd in all_data:
        positions = rd["results"].keys()
        xs, ys = [], []
        for r_au in sorted(positions):
            m = _final_metrics(rd["results"][r_au])
            if np.isnan(m["M_core"]):
                continue
            xs.append(r_au)
            ys.append(m["M_core"])
            ax.scatter(r_au, m["M_core"],
                       color=rd["color"],
                       marker=_MARKERS.get(r_au, "o"),
                       s=_MS**2, zorder=4, edgecolors="white", linewidths=0.5)
        if len(xs) > 1:
            ax.plot(xs, ys, color=rd["color"], lw=1.2, alpha=0.5, zorder=3)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"Posición orbital $a$ [AU]", fontsize=11)
    ax.set_ylabel(r"Masa final del core $M_\mathrm{core}$ [$M_\oplus$]", fontsize=11)
    ax.set_title("Masa final del embrión vs posición orbital — todas las runs",
                 fontsize=13, fontweight="bold", pad=10)

    _add_legend_run_pos(ax, all_data)
    fig.tight_layout()
    _save(fig, savedir, "fig1_mass_radial.pdf")
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Fig 2: f_H2O_final  vs  a [AU]
# ════════════════════════════════════════════════════════════════════════════

def fig_water_radial(all_data: list[dict], savedir: str):
    """Fracción de agua final vs posición orbital."""
    fig, ax = plt.subplots(figsize=(10, 5.5))

    for rd in all_data:
        xs, ys = [], []
        for r_au in sorted(rd["results"].keys()):
            m = _final_metrics(rd["results"][r_au])
            if np.isnan(m["f_H2O"]):
                continue
            xs.append(r_au)
            ys.append(m["f_H2O"])
            ax.scatter(r_au, m["f_H2O"],
                       color=rd["color"],
                       marker=_MARKERS.get(r_au, "o"),
                       s=_MS**2, zorder=4, edgecolors="white", linewidths=0.5)
        if len(xs) > 1:
            ax.plot(xs, ys, color=rd["color"], lw=1.2, alpha=0.5, zorder=3)

    ax.axhline(_WW_THRESHOLD, color="#FFD54F", lw=1.5, ls="--", alpha=0.85,
               label=f"Límite waterworld ({_WW_THRESHOLD:.0f}%)")
    ax.fill_between([min(_POSITIONS_AU)*0.8, max(_POSITIONS_AU)*1.2],
                    _WW_THRESHOLD, ax.get_ylim()[1] if ax.get_ylim()[1] > _WW_THRESHOLD else 100,
                    alpha=0.04, color="#FFD54F", zorder=0)

    ax.set_xscale("log")
    ax.set_xlabel(r"Posición orbital $a$ [AU]", fontsize=11)
    ax.set_ylabel(r"Fracción de $\mathrm{H_2O}$ final [%]", fontsize=11)
    ax.set_title(r"Fracción de agua final vs posición orbital — todas las runs",
                 fontsize=13, fontweight="bold", pad=10)
    ax.set_ylim(bottom=-1)

    _add_legend_run_pos(ax, all_data, extra_handles=[
        Line2D([0], [0], color="#FFD54F", ls="--", lw=1.5,
               label=f"Límite WW ({_WW_THRESHOLD:.0f}%)")
    ])
    fig.tight_layout()
    _save(fig, savedir, "fig2_water_radial.pdf")
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Fig 3 & 4: Grilla multi-panel de evolución temporal
# ════════════════════════════════════════════════════════════════════════════

def fig_temporal_grid(all_data: list[dict], savedir: str,
                      quantity: str = "f_H2O",
                      filename: str = "fig3_fwater_time.pdf"):
    """
    Grilla 2×3 (o 3×2) de paneles.  Un panel por posición orbital.
    Dentro de cada panel: curva temporal para cada run.
    quantity: 'f_H2O' o 'M_core'
    """
    positions = sorted(_POSITIONS_AU)
    npos = len(positions)
    ncols = 3
    nrows = int(np.ceil(npos / ncols))

    fig, axes = plt.subplots(nrows, ncols,
                             figsize=(ncols * 4.5, nrows * 3.4),
                             sharex=False, sharey=False)
    axes = np.array(axes).flatten()

    for idx, r_au in enumerate(positions):
        ax = axes[idx]
        marker = _MARKERS.get(r_au, "o")

        for rd in all_data:
            hist = rd["results"].get(r_au)
            if hist is None or len(hist) == 0:
                continue
            t_yr = hist[:, 0] / YR

            if quantity == "f_H2O":
                M_core = hist[:, 1]
                M_h2o  = hist[:, 2]
                y = np.where(M_core > 1e-30, 100.0 * M_h2o / M_core, 0.0)
                ylabel_panel = r"$f_\mathrm{H_2O}$ [%]"
            else:  # M_core
                y = hist[:, 1] / M_EARTH
                ylabel_panel = r"$M_\mathrm{core}$ [$M_\oplus$]"

            ax.plot(t_yr, y, color=rd["color"], lw=1.6, alpha=0.85)
            # Punto final con marker de posición
            if len(t_yr) > 0:
                ax.scatter(t_yr[-1], y[-1],
                           color=rd["color"], marker=marker,
                           s=50, zorder=5, edgecolors="white", linewidths=0.4)

        if quantity == "f_H2O":
            ax.axhline(_WW_THRESHOLD, color="#FFD54F", lw=1.0, ls="--", alpha=0.7)

        ax.set_xscale("log")
        if quantity == "M_core":
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                ax.set_yscale("log")
        ax.set_title(f"a = {r_au:.1f} AU  [{marker}]", fontsize=10, pad=4)
        ax.set_xlabel("Tiempo [años]", fontsize=8)
        if idx % ncols == 0:
            ax.set_ylabel(ylabel_panel, fontsize=9)

    # Ocultar paneles sobrantes
    for extra in axes[npos:]:
        extra.set_visible(False)

    # Leyenda de runs en la última celda disponible
    if npos < len(axes):
        ax_leg = axes[npos]
        ax_leg.set_visible(True)
        ax_leg.axis("off")
        handles = [Line2D([0], [0], color=rd["color"], lw=2.0,
                          label=_pretty(rd["label"])) for rd in all_data]
        ax_leg.legend(handles=handles, loc="center", fontsize=8,
                      framealpha=0.3, title="Runs", title_fontsize=9)
    else:
        # Leyenda en panel superior derecho
        handles = [Line2D([0], [0], color=rd["color"], lw=2.0,
                          label=_pretty(rd["label"])) for rd in all_data]
        fig.legend(handles=handles, loc="upper right", fontsize=7,
                   framealpha=0.3, title="Runs", ncols=2)

    suptitle = (r"Fracción de $\mathrm{H_2O}$ — evolución temporal por posición orbital"
                if quantity == "f_H2O"
                else r"Masa del core — evolución temporal por posición orbital")
    fig.suptitle(suptitle, fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    _save(fig, savedir, filename)
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Fig 5: Stacked-bar agrupado por posición orbital
# ════════════════════════════════════════════════════════════════════════════

def fig_architecture_grouped(all_data: list[dict], savedir: str):
    """
    Barras agrupadas: grupos = posiciones orbitales, dentro de cada grupo
    una barra por run (coloreada).  Cada barra apilada Silicatos / CO2 / H2O.
    """
    positions = sorted(_POSITIONS_AU)
    npos   = len(positions)
    nruns  = len(all_data)

    group_width = 0.8
    bar_width   = group_width / nruns
    fig_w = max(14, npos * nruns * 0.28 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5.5))

    for gi, r_au in enumerate(positions):
        for ri, rd in enumerate(all_data):
            hist  = rd["results"].get(r_au)
            m     = _final_metrics(hist)
            if np.isnan(m["M_core"]):
                continue

            x_center = gi + (ri - nruns / 2 + 0.5) * bar_width
            M_tot    = m["M_core"]
            M_h2o    = hist[-1, 2] / M_EARTH
            M_co2    = hist[-1, 3] / M_EARTH
            M_sil    = max(0.0, M_tot - M_h2o - M_co2)

            ax.bar(x_center, M_sil, width=bar_width * 0.88, bottom=0,
                   color=_COMP_COLS["silicates"], edgecolor="black", linewidth=0.2, zorder=3)
            ax.bar(x_center, M_co2, width=bar_width * 0.88, bottom=M_sil,
                   color=_COMP_COLS["CO2"],       edgecolor="black", linewidth=0.2, zorder=3)
            ax.bar(x_center, M_h2o, width=bar_width * 0.88, bottom=M_sil + M_co2,
                   color=_COMP_COLS["H2O"],       edgecolor="black", linewidth=0.2, zorder=3)

            # Marcar waterworld con borde de color de run
            if m["f_H2O"] > _WW_THRESHOLD:
                ax.bar(x_center, M_tot, width=bar_width * 0.88,
                       fill=False, edgecolor=rd["color"], linewidth=1.8, zorder=5)

    # Ticks grupos (posición central de cada grupo)
    ax.set_xticks(np.arange(npos))
    ax.set_xticklabels([f"{r:.1f} AU" for r in positions], fontsize=9)
    ax.set_ylabel(r"Masa final [$M_\oplus$]", fontsize=11)
    ax.set_title("Arquitectura final del core agrupado por posición orbital",
                 fontsize=13, fontweight="bold", pad=10)

    # Leyenda compuesta: composición + runs
    comp_handles = [Patch(color=_COMP_COLS["silicates"], label="Silicatos"),
                    Patch(color=_COMP_COLS["CO2"],       label=r"$\mathrm{CO_2}$"),
                    Patch(color=_COMP_COLS["H2O"],       label=r"$\mathrm{H_2O}$")]
    run_handles  = [Patch(facecolor="none", edgecolor=rd["color"],
                          linewidth=1.8, label=_pretty(rd["label"])) for rd in all_data]

    leg1 = ax.legend(handles=comp_handles, loc="upper left",  fontsize=8, title="Composición")
    ax.legend(handles=run_handles, loc="upper right",
              fontsize=7, title="Runs (borde = WW)", ncols=2)
    ax.add_artist(leg1)

    fig.tight_layout()
    _save(fig, savedir, "fig5_architecture_grouped.pdf")
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Fig 6: Heatmap 2D  runs × posiciones
# ════════════════════════════════════════════════════════════════════════════

def fig_heatmap(all_data: list[dict], savedir: str):
    """
    Heatmap 2D: filas = runs, columnas = posiciones orbitales.
    Color = f_H2O_final (%).  Texto en cada celda = M_final [M⊕].
    """
    positions = sorted(_POSITIONS_AU)
    nruns     = len(all_data)
    npos      = len(positions)

    matrix_fh2o = np.full((nruns, npos), np.nan)
    matrix_mass  = np.full((nruns, npos), np.nan)

    for ri, rd in enumerate(all_data):
        for ci, r_au in enumerate(positions):
            hist = rd["results"].get(r_au)
            m    = _final_metrics(hist)
            matrix_fh2o[ri, ci] = m["f_H2O"]
            matrix_mass[ri,  ci] = m["M_core"]

    fig, ax = plt.subplots(figsize=(max(9, npos * 1.5), max(5, nruns * 0.8)))

    cmap  = matplotlib.colormaps.get_cmap("YlOrRd")
    cmap.set_bad("#2A2A2A")
    vmax  = max(_WW_THRESHOLD * 3, np.nanmax(matrix_fh2o) * 1.05)

    im = ax.imshow(matrix_fh2o, aspect="auto", cmap=cmap,
                   vmin=0, vmax=vmax, interpolation="nearest")

    # Texto en cada celda
    for ri in range(nruns):
        for ci in range(npos):
            fh2o = matrix_fh2o[ri, ci]
            mass = matrix_mass[ri, ci]
            if np.isnan(fh2o):
                cell_text = "N/A"
                tc = "white"
            else:
                cell_text = f"{mass:.2f}\n{fh2o:.1f}%"
                # Contraste automático
                tc = "black" if fh2o < vmax * 0.60 else "white"
            ax.text(ci, ri, cell_text, ha="center", va="center",
                    fontsize=7.5, color=tc, fontweight="bold")

    # Línea de umbral waterworld (en las celdas > 10%)
    for ri in range(nruns):
        for ci in range(npos):
            if not np.isnan(matrix_fh2o[ri, ci]) and matrix_fh2o[ri, ci] >= _WW_THRESHOLD:
                rect = plt.Rectangle((ci - 0.5, ri - 0.5), 1, 1,
                                     fill=False, edgecolor="#FFD54F",
                                     linewidth=2.0, zorder=5)
                ax.add_patch(rect)

    # Ejes
    ax.set_xticks(np.arange(npos))
    ax.set_xticklabels([f"{r:.1f} AU" for r in positions], fontsize=9)
    ax.set_yticks(np.arange(nruns))
    ax.set_yticklabels([_pretty(rd["label"]) for rd in all_data], fontsize=8)
    ax.set_xlabel("Posición orbital del embrión [AU]", fontsize=11)
    ax.set_title(
        r"Heatmap: $f_\mathrm{H_2O}$ final [%] — todas las runs × todas las posiciones"
        + f"\n(Marco dorado = Waterworld, >{_WW_THRESHOLD:.0f}%  |  texto: masa [$M_\\oplus$]%)",
        fontsize=12, fontweight="bold", pad=10,
    )

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(r"$f_\mathrm{H_2O}$ final [%]", fontsize=10)
    cbar.ax.axhline(_WW_THRESHOLD, color="#FFD54F", lw=1.5, ls="--")

    fig.tight_layout()
    _save(fig, savedir, "fig6_heatmap_runs_positions.pdf")
    return fig


# ════════════════════════════════════════════════════════════════════════════
# Helpers de leyenda
# ════════════════════════════════════════════════════════════════════════════

def _add_legend_run_pos(ax, all_data: list[dict],
                        extra_handles: list | None = None):
    """
    Leyenda compuesta:
      Columna 1: colores de run
      Columna 2: marcadores de posición
    """
    run_handles = [Line2D([0], [0], color=rd["color"], marker="o",
                          markersize=7, lw=1.5,
                          label=_pretty(rd["label"])) for rd in all_data]
    pos_handles = [Line2D([0], [0], color="gray",
                          marker=_MARKERS[r], markersize=8, lw=0,
                          label=f"{r:.1f} AU") for r in _POSITIONS_AU]

    all_handles = run_handles + [
        Line2D([0], [0], lw=0)  # separador invisible
    ] + pos_handles

    if extra_handles:
        all_handles = extra_handles + all_handles

    ax.legend(handles=all_handles, fontsize=7.5,
              ncols=2 if len(run_handles) > 5 else 1,
              loc="best", framealpha=0.3,
              handlelength=1.8, handletextpad=0.6)


# ════════════════════════════════════════════════════════════════════════════
# Pipeline principal
# ════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    dataroot:  str   = "data/gaps_pipeline/t_5e6",
    savedir:   str   = "figures/pa3/multi_embryo",
    t_min_yr:  float = 100.0,
    positions: list  | None = None,
    show:      bool  = True,
):
    """
    Punto de entrada principal.

    Para cada run en `dataroot/`, ejecuta PA3 de forma INDEPENDIENTE para
    cada embrión en `positions` y genera 6 figuras comparativas.

    Parameters
    ----------
    dataroot  : directorio raíz con subdirectorios de runs
    savedir   : directorio de salida (se crea si no existe)
    t_min_yr  : tiempo mínimo de snapshots a usar [años]
    positions : lista de posiciones orbitales [AU]  (default: global _POSITIONS_AU)
    show      : si True llama plt.show() al final
    """
    if positions is None:
        positions = _POSITIONS_AU

    print("=" * 72)
    print(f"  Multi-Embryo Pipeline — dataroot : {dataroot}")
    print(f"  Posiciones [AU] : {positions}")
    print(f"  t_min_yr        : {t_min_yr}")
    print(f"  savedir         : {savedir}")
    print("=" * 72)

    runs = _collect_runs(dataroot)
    if not runs:
        print(f"[ERROR] No se encontraron runs con HDF5 en {dataroot}")
        return

    print(f"\n  Runs encontradas ({len(runs)}):")
    for name, _ in runs:
        print(f"    • {name}")
    print()

    all_data = []
    for idx, (name, path) in enumerate(runs):
        color = _PALETTE[idx % len(_PALETTE)]
        print(f"\n{'─' * 64}")
        print(f"  [{idx+1}/{len(runs)}] {name}")
        print(f"{'─' * 64}")
        try:
            results = _run_independent(path, positions, t_min_yr)
            all_data.append({
                "label":   name,
                "color":   color,
                "results": results,
            })
            # Mini resumen
            for r_au in sorted(positions):
                m = _final_metrics(results.get(r_au))
                tipo = "💧 WW" if not np.isnan(m["f_H2O"]) and m["f_H2O"] > _WW_THRESHOLD else "🪨  Roc"
                mstr = f"{m['M_core']:.3f}" if not np.isnan(m["M_core"]) else " N/A "
                fstr = f"{m['f_H2O']:.1f}" if not np.isnan(m["f_H2O"]) else " N/A "
                print(f"    a={r_au:4.1f} AU → M={mstr} M⊕   f_H2O={fstr}%  {tipo}")
        except Exception as exc:
            print(f"  [ERROR] {name}: {exc}")
            import traceback; traceback.print_exc()
            continue

    if not all_data:
        print("\n[ERROR] Ninguna run fue procesada.")
        return

    print(f"\n{'=' * 72}")
    print("  Generando figuras…")
    print(f"{'=' * 72}\n")

    fig_mass_radial(all_data, savedir)
    fig_water_radial(all_data, savedir)
    fig_temporal_grid(all_data, savedir, quantity="f_H2O",  filename="fig3_fwater_time.pdf")
    fig_temporal_grid(all_data, savedir, quantity="M_core", filename="fig4_mass_time.pdf")
    fig_architecture_grouped(all_data, savedir)
    fig_heatmap(all_data, savedir)

    print(f"\n  ✓ Figuras guardadas en:  {savedir}/")
    for fn in ["fig1_mass_radial.pdf", "fig2_water_radial.pdf",
               "fig3_fwater_time.pdf", "fig4_mass_time.pdf",
               "fig5_architecture_grouped.pdf", "fig6_heatmap_runs_positions.pdf"]:
        print(f"    • {fn}")

    if show:
        plt.show()


# ════════════════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Pipeline multi-embrión PA3 — todas las runs t_5e6."
    )
    parser.add_argument("--dataroot", default="data/gaps_pipeline/t_5e6",
                        help="Directorio raíz con las carpetas de runs")
    parser.add_argument("--savedir",  default="figures/pa3/multi_embryo",
                        help="Directorio de salida para los PDFs")
    parser.add_argument("--t_min_yr", type=float, default=100.0,
                        help="Tiempo mínimo de snapshot [años]")
    parser.add_argument("--positions", type=float, nargs="+",
                        default=None,
                        help="Posiciones orbitales [AU] (default: 1 2 2.5 3 5 10)")
    parser.add_argument("--no-show", action="store_true",
                        help="No llamar plt.show()")
    args = parser.parse_args()

    run_pipeline(
        dataroot  = args.dataroot,
        savedir   = args.savedir,
        t_min_yr  = args.t_min_yr,
        positions = args.positions,
        show      = not args.no_show,
    )

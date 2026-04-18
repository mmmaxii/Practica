"""
plot_pebble_runs.py
===================
Pipeline de plots comparativos para las runs de pebble accretion en t_5e6.

Para cada subdirectorio en `data_gaps_pipeline/t_5e6/`:
  1. Carga los snapshots HDF5 con PebbleAccretion3
  2. Coloca un embrión en 2 AU y corre el modelo de acreción
  3. Superpone todas las runs en figuras comparativas:

     Figura 1 — Masa del core M(t) para todas las runs
     Figura 2 — Fracción de agua f_H2O(t) para todas las runs
     Figura 3 — Arquitectura del core al final: stacked-bar chart
                  silicatos / CO2 / H2O  por run

Uso:
    python plot_pebble_runs.py [--dataroot data_gaps_pipeline/t_5e6]
                               [--savedir figs_pa3]
                               [--r_planet 2.0]
                               [--t_min_yr 100]
"""

import os
import sys
import glob
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import dustpy.constants as c

# ── Importar el módulo de acreción ──────────────────────────────────────────
# Añadimos el directorio padre de PA3Py al path si es necesario
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3

# ── Estética global ──────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.titlesize":  13,
    "axes.labelsize":  11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.grid":       True,
    "grid.alpha":      0.15,
    "figure.dpi":      130,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
})

# Paleta de colores — distinguibles para hasta ~10 runs
_PALETTE = [
    "#4FC3F7",  # azul claro
    "#FF8A65",  # naranja salmón
    "#81C784",  # verde claro
    "#CE93D8",  # lila
    "#FFCC80",  # amarillo cálido
    "#80CBC4",  # teal
    "#F48FB1",  # rosa
    "#AED581",  # lima
    "#90CAF9",  # azul pálido
    "#FFAB91",  # melocotón
]

# Colores de composición del core (stacked bar)
_COMP_COLS = {
    "silicates": "#EF9A9A",
    "CO2":       "#80CBC4",
    "H2O":       "#4FC3F7",
}
_COMP_LABELS = {
    "silicates": "Silicatos",
    "CO2":       r"$\mathrm{CO_2}$",
    "H2O":       r"$\mathrm{H_2O}$",
}

M_EARTH = c.M_earth   # g
AU      = c.au        # cm
YR      = c.year      # s


# ════════════════════════════════════════════════════════════════════════════
# Utilidades
# ════════════════════════════════════════════════════════════════════════════

def _pretty_label(run_name: str) -> str:
    """Convierte nombre de carpeta en label legible para la leyenda."""
    return run_name.replace("_", " ")


def _collect_runs(dataroot: str) -> list[tuple[str, str]]:
    """
    Devuelve lista de (run_name, run_path) ordenada alfabéticamente.
    Solo incluye subdirectorios que contengan al menos un archivo data*.hdf5.
    """
    runs = []
    for entry in sorted(os.listdir(dataroot)):
        path = os.path.join(dataroot, entry)
        if not os.path.isdir(path):
            continue
        hdf5_files = glob.glob(os.path.join(path, "data*.hdf5"))
        if hdf5_files:
            runs.append((entry, path))
        else:
            print(f"  [skip] {entry}  (sin archivos data*.hdf5)")
    return runs


def _run_pa3(run_path: str, r_planet_au: float, t_min_yr: float):
    """
    Carga el run y ejecuta PA3 para un planeta en r_planet_au.

    Returns
    -------
    pa3   : PebbleAccretionModule3
    hist  : np.ndarray  shape (Nt_acret, 7)
               cols: [t_s, M_core_g, M_H2O_g, M_CO2_g, M_sil_g, r_snow_AU, M_iso_g]
    """
    pa3  = PebbleAccretionModule3.from_datadir(run_path, t_min_yr=t_min_yr)
    results = pa3.run_growth([r_planet_au])
    hist = results[r_planet_au]   # (Nt_acret, 7)
    return pa3, hist


# ════════════════════════════════════════════════════════════════════════════
# Figura 1: Masa del core M(t) — todas las runs
# ════════════════════════════════════════════════════════════════════════════

def plot_mass_evolution(
    run_data: list[dict],
    savedir:  str | None = None,
    r_planet_au: float   = 2.0,
) -> tuple:
    """
    Líneas M_core(t) [M_⊕] para todas las runs en un mismo panel.

    `run_data` : lista de dicts con keys 'label', 'hist', 'color'
    """
    fig, ax = plt.subplots(figsize=(10, 5.5))

    for rd in run_data:
        hist  = rd["hist"]
        if len(hist) == 0:
            print(f"  [WARN] {rd['label']}: sin historial de acreción.")
            continue
        t_yr   = hist[:, 0] / YR
        M_core = hist[:, 1] / M_EARTH   # M_⊕
        M_iso  = hist[:, 6] / M_EARTH   # M_iso_peb en M_⊕

        ax.plot(t_yr, M_core, color=rd["color"], lw=2.0,
                label=_pretty_label(rd["label"]))
        # Línea punteada de M_iso (misma run, misma línea pero alfa reducida)
        ax.plot(t_yr, M_iso, color=rd["color"], lw=1.0, ls=":",
                alpha=0.45)

    # Línea fantasma para la leyenda de M_iso
    ax.plot([], [], "gray", lw=1.0, ls=":", alpha=0.6,
            label=r"$M_\mathrm{iso}$ (punteada, por run)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Tiempo [años]", fontsize=11)
    ax.set_ylabel(r"Masa del core [$M_\oplus$]", fontsize=11)
    ax.set_title(
        rf"Evolución de masa del embrión @ {r_planet_au:.1f} AU — todas las runs",
        fontsize=13, fontweight="bold", pad=10,
    )

    # Anotación de masa de aislamiento
    ax.legend(fontsize=8, framealpha=0.3, loc="upper left")
    fig.tight_layout()

    if savedir:
        os.makedirs(savedir, exist_ok=True)
        out = os.path.join(savedir, "pebble_mass_evolution.pdf")
        fig.savefig(out, bbox_inches="tight")
        print(f"  → Guardado: {out}")

    return fig, ax


# ════════════════════════════════════════════════════════════════════════════
# Figura 2: Fracción de agua f_H2O(t) — todas las runs
# ════════════════════════════════════════════════════════════════════════════

def plot_water_fraction(
    run_data:    list[dict],
    savedir:     str | None = None,
    r_planet_au: float      = 2.0,
    threshold_ww: float     = 0.10,   # fracción > 10 % → waterworld
) -> tuple:
    """
    Líneas f_H2O(t) [%] para todas las runs.
    Línea horizontal punteada en threshold_ww*100 como referencia waterworld.
    """
    fig, ax = plt.subplots(figsize=(10, 5.5))

    for rd in run_data:
        hist = rd["hist"]
        if len(hist) == 0:
            continue
        t_yr  = hist[:, 0] / YR
        M_core = hist[:, 1]               # g
        M_h2o  = hist[:, 2]               # g  (col 2 = H2O)
        # Fracción en porcentaje; evita división por cero
        f_h2o  = np.where(M_core > 1e-30, 100.0 * M_h2o / M_core, 0.0)

        ax.plot(t_yr, f_h2o, color=rd["color"], lw=2.0,
                label=_pretty_label(rd["label"]))

    # Línea referencia waterworld
    ax.axhline(threshold_ww * 100,
               color="white", lw=1.2, ls="--", alpha=0.6,
               label=rf"Límite waterworld ({threshold_ww*100:.0f}%)")

    ax.set_xscale("log")
    ax.set_xlabel("Tiempo [años]", fontsize=11)
    ax.set_ylabel(r"Fracción de $\mathrm{H_2O}$ en el core [%]", fontsize=11)
    ax.set_title(
        rf"Evolución de fracción de agua @ {r_planet_au:.1f} AU — todas las runs",
        fontsize=13, fontweight="bold", pad=10,
    )
    ax.set_ylim(bottom=-1)

    ax.legend(fontsize=8, framealpha=0.3, loc="upper left")
    fig.tight_layout()

    if savedir:
        os.makedirs(savedir, exist_ok=True)
        out = os.path.join(savedir, "pebble_water_fraction.pdf")
        fig.savefig(out, bbox_inches="tight")
        print(f"  → Guardado: {out}")

    return fig, ax


# ════════════════════════════════════════════════════════════════════════════
# Figura 3: Arquitectura del core — stacked-bar chart (estado FINAL)
# ════════════════════════════════════════════════════════════════════════════

def plot_core_architecture(
    run_data:    list[dict],
    savedir:     str | None = None,
    r_planet_au: float      = 2.0,
) -> tuple:
    """
    Gráfico de barras horizontales apiladas con la composición final del core:
        [silicatos | CO2 | H2O]  en M_⊕

    Cada barra es una run.  Las masas están en escala lineal para visualizar
    la arquitectura interna real.
    """
    labels = []
    M_sil_f, M_co2_f, M_h2o_f, M_tot_f = [], [], [], []

    for rd in run_data:
        hist = rd["hist"]
        if len(hist) == 0:
            continue
        row = hist[-1]          # último snapshot
        M_tot = row[1]          # g
        M_h2o = row[2]          # g
        M_co2 = row[3]          # g
        M_sil = max(0.0, M_tot - M_h2o - M_co2)  # silicatos  = resto

        labels.append(_pretty_label(rd["label"]))
        M_tot_f.append(M_tot / M_EARTH)
        M_h2o_f.append(M_h2o / M_EARTH)
        M_co2_f.append(M_co2 / M_EARTH)
        M_sil_f.append(M_sil / M_EARTH)

    n = len(labels)
    if n == 0:
        print("  [WARN] No hay runs con historial para plot_core_architecture.")
        return None, None

    y_pos = np.arange(n)
    bar_h = 0.55

    fig, ax = plt.subplots(figsize=(11, max(4, 0.9 * n + 2)))

    bars_sil = ax.barh(y_pos, M_sil_f, height=bar_h, color=_COMP_COLS["silicates"],
                       label=_COMP_LABELS["silicates"], zorder=3)
    bars_co2 = ax.barh(y_pos, M_co2_f, height=bar_h, left=M_sil_f,
                       color=_COMP_COLS["CO2"], label=_COMP_LABELS["CO2"], zorder=3)
    left_h2o = np.array(M_sil_f) + np.array(M_co2_f)
    bars_h2o = ax.barh(y_pos, M_h2o_f, height=bar_h, left=left_h2o,
                       color=_COMP_COLS["H2O"], label=_COMP_LABELS["H2O"], zorder=3)

    # Anotar masa total al final de cada barra
    for i, (m_tot, m_sil, m_co2, m_h2o) in \
            enumerate(zip(M_tot_f, M_sil_f, M_co2_f, M_h2o_f)):
        f_h2o = 100 * m_h2o / (m_tot + 1e-30)
        tipo  = "💧 WW" if f_h2o > 10 else "🪨  Roc"
        ax.text(m_tot + 0.03, y_pos[i],
                f"{m_tot:.2f} $M_\\oplus$  {tipo}",
                va="center", ha="left", fontsize=8.5,
                color="white", zorder=4)

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r"Masa [$M_\oplus$]", fontsize=11)
    ax.set_title(
        rf"Arquitectura del core @ {r_planet_au:.1f} AU — estado final",
        fontsize=13, fontweight="bold", pad=10,
    )
    ax.legend(loc="lower right", fontsize=9, framealpha=0.25)
    ax.invert_yaxis()   # primer run arriba
    ax.set_xlim(left=0)
    ax.xaxis.set_minor_locator(mticker.AutoMinorLocator())

    fig.tight_layout()

    if savedir:
        os.makedirs(savedir, exist_ok=True)
        out = os.path.join(savedir, "pebble_core_architecture.pdf")
        fig.savefig(out, bbox_inches="tight")
        print(f"  → Guardado: {out}")

    return fig, ax


# ════════════════════════════════════════════════════════════════════════════
# Figura 4 (bonus): Panel combinado M(t) + f_H2O(t), 2 filas × 1 col
# ════════════════════════════════════════════════════════════════════════════

def plot_combined(
    run_data:    list[dict],
    savedir:     str | None = None,
    r_planet_au: float      = 2.0,
    threshold_ww: float     = 0.10,
) -> tuple:
    """
    Panel combinado de masa (arriba) y fracción de agua (abajo), eje X compartido.
    """
    fig, (ax_m, ax_w) = plt.subplots(2, 1, figsize=(10, 9),
                                      sharex=True, gridspec_kw={"hspace": 0.08})

    for rd in run_data:
        hist = rd["hist"]
        if len(hist) == 0:
            continue
        t_yr   = hist[:, 0] / YR
        M_core = hist[:, 1] / M_EARTH
        M_h2o  = hist[:, 2]
        M_iso  = hist[:, 6] / M_EARTH
        f_h2o  = np.where(hist[:, 1] > 1e-30, 100.0 * M_h2o / hist[:, 1], 0.0)

        ax_m.plot(t_yr, M_core, color=rd["color"], lw=2.0,
                  label=_pretty_label(rd["label"]))
        ax_m.plot(t_yr, M_iso, color=rd["color"], lw=0.9, ls=":", alpha=0.40)

        ax_w.plot(t_yr, f_h2o, color=rd["color"], lw=2.0)

    # Leyenda solo en panel superior (para no duplicar)
    ax_m.legend(fontsize=8, framealpha=0.25, loc="upper left",
                ncols=min(3, len(run_data)))
    ax_m.plot([], [], "gray", lw=0.9, ls=":", alpha=0.6,
              label=r"$M_\mathrm{iso}$ (punteada)")
    ax_m.set_xscale("log")
    ax_m.set_yscale("log")
    ax_m.set_ylabel(r"$M_\mathrm{core}$ [$M_\oplus$]", fontsize=11)
    ax_m.set_title(
        rf"Pebble Accretion @ {r_planet_au:.1f} AU — todas las runs",
        fontsize=13, fontweight="bold", pad=10,
    )

    ax_w.axhline(threshold_ww * 100, color="white",
                 lw=1.2, ls="--", alpha=0.55,
                 label=rf"WW límite ({threshold_ww*100:.0f}%)")
    ax_w.legend(fontsize=8, framealpha=0.25, loc="upper left")
    ax_w.set_xlabel("Tiempo [años]", fontsize=11)
    ax_w.set_ylabel(r"$f_{\mathrm{H_2O}}$ [%]", fontsize=11)
    ax_w.set_ylim(bottom=-1)

    fig.tight_layout()

    if savedir:
        os.makedirs(savedir, exist_ok=True)
        out = os.path.join(savedir, "pebble_combined.pdf")
        fig.savefig(out, bbox_inches="tight")
        print(f"  → Guardado: {out}")

    return fig, (ax_m, ax_w)


# ════════════════════════════════════════════════════════════════════════════
# Pipeline principal
# ════════════════════════════════════════════════════════════════════════════

def run_pipeline(
    dataroot:    str   = "data_gaps_pipeline/t_5e6",
    savedir:     str   = "figs_pa3",
    r_planet_au: float = 2.0,
    t_min_yr:    float = 100.0,
    show:        bool  = True,
):
    """
    Punto de entrada principal.

    1. Descubre todas las runs con HDF5 en `dataroot/`.
    2. Para cada run, ejecuta PA3 con un embrión en `r_planet_au` AU.
    3. Genera las 4 figuras comparativas y las guarda en `savedir/`.

    Parameters
    ----------
    dataroot    : directorio raíz con subdirectorios de runs
    savedir     : directorio de salida para los PDFs
    r_planet_au : posición del embrión en AU
    t_min_yr    : tiempo mínimo (años) a usar en el modelo
    show        : si True llama plt.show() al final
    """
    print("=" * 70)
    print(f"  Pipeline PA3 — runs en: {dataroot}")
    print(f"  Planeta @ {r_planet_au} AU  |  t_min = {t_min_yr:.0f} yr")
    print("=" * 70)

    runs = _collect_runs(dataroot)
    if not runs:
        print(f"[ERROR] No se encontraron runs con HDF5 en {dataroot}")
        return

    print(f"\n  Runs encontradas ({len(runs)}):")
    for name, path in runs:
        print(f"    • {name}")
    print()

    run_data = []
    for idx, (name, path) in enumerate(runs):
        color = _PALETTE[idx % len(_PALETTE)]
        print(f"\n{'─'*60}")
        print(f"  [{idx+1}/{len(runs)}] Procesando: {name}")
        print(f"{'─'*60}")
        try:
            pa3, hist = _run_pa3(path, r_planet_au, t_min_yr)
            pa3.summary({r_planet_au: hist})
            run_data.append({
                "label": name,
                "hist":  hist,
                "pa3":   pa3,
                "color": color,
            })
        except Exception as e:
            print(f"  [ERROR] {name}: {e}")
            import traceback; traceback.print_exc()
            continue

    if not run_data:
        print("[ERROR] Ninguna run fue procesada correctamente.")
        return

    print(f"\n{'='*70}")
    print("  Generando figuras...")
    print(f"{'='*70}\n")

    # Figuras individuales
    plot_mass_evolution(run_data, savedir=savedir, r_planet_au=r_planet_au)
    plot_water_fraction(run_data, savedir=savedir, r_planet_au=r_planet_au)
    plot_core_architecture(run_data, savedir=savedir, r_planet_au=r_planet_au)

    # Panel combinado (bonus)
    plot_combined(run_data, savedir=savedir, r_planet_au=r_planet_au)

    print(f"\n  ✓ Figuras guardadas en: {savedir}/")
    print(f"    • pebble_mass_evolution.pdf")
    print(f"    • pebble_water_fraction.pdf")
    print(f"    • pebble_core_architecture.pdf")
    print(f"    • pebble_combined.pdf")

    if show:
        plt.show()


# ════════════════════════════════════════════════════════════════════════════
# CLI
# ════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Pipeline de plots PA3 comparativos para runs de t_5e6."
    )
    parser.add_argument(
        "--dataroot", default="data_gaps_pipeline/t_5e6",
        help="Directorio raíz con las carpetas de runs (default: data_gaps_pipeline/t_5e6)"
    )
    parser.add_argument(
        "--savedir", default="figs_pa3",
        help="Directorio de salida para los PDFs (default: figs_pa3)"
    )
    parser.add_argument(
        "--r_planet", type=float, default=2.0,
        help="Posición del embrión en AU (default: 2.0)"
    )
    parser.add_argument(
        "--t_min_yr", type=float, default=100.0,
        help="Tiempo mínimo en años para usar snapshots (default: 100)"
    )
    parser.add_argument(
        "--no-show", action="store_true",
        help="No llamar plt.show() (útil en scripts no-interactivos)"
    )
    args = parser.parse_args()

    run_pipeline(
        dataroot    = args.dataroot,
        savedir     = args.savedir,
        r_planet_au = args.r_planet,
        t_min_yr    = args.t_min_yr,
        show        = not args.no_show,
    )

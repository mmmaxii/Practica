"""
run_diagnostics.py — Script de post-procesamiento multi-run
=============================================================

Define múltiples configuraciones de simulación (sin gap, Kanagawa, Duffell),
las corre secuencialmente con su propio datadir, y genera figuras diagnóstico
comparativas usando SnowlineDiagnostics.

Configuraciones disponibles
----------------------------
  "baseline"         → sin gap planetario
  "kanagawa_saturn"  → gap Kanagawa, Saturno (~95 M⊕) en 9.5 AU
  "kanagawa_jupiter" → gap Kanagawa, Júpiter (~318 M⊕) en 5.2 AU
  "duffell_neptune"  → gap Duffell, Neptuno (~17 M⊕) en 30 AU
  "duffell_jupiter"  → gap Duffell, Júpiter en 5.2 AU con imprint en Σ

Uso
---
    python run_diagnostics.py            → corre todos los runs + diagnósticos
    python run_diagnostics.py --diag     → solo diagnósticos (no simula de nuevo)
    python run_diagnostics.py --compare  → solo figuras de comparación

Requiere
--------
    pipeline_snowlines.py  (WaterworldPipeline)
    plot_diagnostics.py    (SnowlineDiagnostics)
    pressure_bumps.py      (PressureBumpsMixin — ya integrado en el pipeline)
"""

import os
import sys
import time
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import dustpy.constants as c
from dustpy import hdf5writer

from pipeline_snowlines import WaterworldPipeline
from plot_diagnostics    import SnowlineDiagnostics

# ══════════════════════════════════════════════════════════════════════════════
# Directorio raíz de salida
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR    = "data_gaps_pipeline"
FIGURES_DIR = "figs_gaps_pipeline"

# ══════════════════════════════════════════════════════════════════════════════
# Configuraciones de cada run
# ══════════════════════════════════════════════════════════════════════════════
#
# Cada entrada es un dict con:
#   label        : nombre corto para identificar el run
#   datadir      : directorio HDF5 del run
#   gap_model    : "none" | "kanagawa" | "duffell"
#   M_planet     : masa del planeta [g — CGS] (ignorado si gap_model="none")
#   a_planet_au  : semieje del planeta [AU]
#   alpha_ref    : α de referencia [adim] (default: 1e-3)
#   imprint      : bool — aplicar gap a Σ desde t=0
#   active_species: lista de volátiles
#   t_end_years  : tiempo final de integración [yr]
#   num_snapshots: número de snapshots
#   M_star_Msun  : masa estelar [M☉]
#   Nr           : celdas radiales
#   rmax_au      : radio exterior [AU]
#
RUNS = [
    {
        "label":          "baseline",
        "datadir":        os.path.join(BASE_DIR, "baseline"),
        "gap_model":      "none",
        "M_planet":       0.0,
        "a_planet_au":    0.0,
        "alpha_ref":      1e-3,
        "imprint":        False,
        "active_species": ["H2O"],
        "t_end_years":    1e6,
        "num_snapshots":  50,
        "M_star_Msun":    1.0,
        "Nr":             200,
        "rmax_au":        300.0,
    },
    {
        "label":          "kanagawa_saturn",
        "datadir":        os.path.join(BASE_DIR, "kanagawa_saturn"),
        "gap_model":      "kanagawa",
        "M_planet":       95.159 * c.M_earth,   # Saturno
        "a_planet_au":    9.5826,
        "alpha_ref":      1e-3,
        "imprint":        False,
        "active_species": ["H2O"],
        "t_end_years":    1e6,
        "num_snapshots":  50,
        "M_star_Msun":    1.0,
        "Nr":             200,
        "rmax_au":        300.0,
    },
    {
        "label":          "kanagawa_jupiter",
        "datadir":        os.path.join(BASE_DIR, "kanagawa_jupiter"),
        "gap_model":      "kanagawa",
        "M_planet":       317.8 * c.M_earth,    # Júpiter
        "a_planet_au":    5.2038,
        "alpha_ref":      1e-3,
        "imprint":        True,                  # gap impreso en Σ desde t=0
        "active_species": ["H2O"],
        "t_end_years":    1e6,
        "num_snapshots":  50,
        "M_star_Msun":    1.0,
        "Nr":             200,
        "rmax_au":        300.0,
    },
    {
        "label":          "duffell_neptune",
        "datadir":        os.path.join(BASE_DIR, "duffell_neptune"),
        "gap_model":      "duffell",
        "M_planet":       17.147 * c.M_earth,   # Neptuno
        "a_planet_au":    30.07,
        "alpha_ref":      1e-3,
        "imprint":        False,
        "active_species": ["H2O"],
        "t_end_years":    1e6,
        "num_snapshots":  50,
        "M_star_Msun":    1.0,
        "Nr":             200,
        "rmax_au":        300.0,
    },
    {
        "label":          "duffell_jupiter_imprint",
        "datadir":        os.path.join(BASE_DIR, "duffell_jupiter_imprint"),
        "gap_model":      "duffell",
        "M_planet":       317.8 * c.M_earth,    # Júpiter
        "a_planet_au":    5.2038,
        "alpha_ref":      1e-3,
        "imprint":        True,
        "active_species": ["H2O"],
        "t_end_years":    1e6,
        "num_snapshots":  50,
        "M_star_Msun":    1.0,
        "Nr":             200,
        "rmax_au":        300.0,
    },
]


# ══════════════════════════════════════════════════════════════════════════════
# Función: construir y correr un run
# ══════════════════════════════════════════════════════════════════════════════

def build_and_run(cfg: dict):
    """
    Construye y corre un WaterworldPipeline con la configuración dada.

    Parameters
    ----------
    cfg : dict
        Configuración del run (ver RUNS arriba).

    Returns
    -------
    datadir : str
        Directorio donde se guardaron los HDF5.
    elapsed : float
        Tiempo de ejecución en segundos.
    """
    label   = cfg["label"]
    datadir = cfg["datadir"]

    sep = "═" * 60
    print(f"\n{sep}")
    print(f"  RUN: {label}")
    print(f"  datadir: {datadir}")
    print(sep)

    t0 = time.time()

    pipeline = WaterworldPipeline(datadir=datadir)

    # ── 1. Configuración de especie y grilla ──────────────────────────────────
    pipeline.active_species = cfg["active_species"]
    pipeline.setup_grid(
        rmin = 1.0 * c.au,
        rmax = cfg["rmax_au"] * c.au,
        Nr   = cfg["Nr"],
    )
    pipeline.setup_star(M_star_Msun=cfg["M_star_Msun"])

    # ── 2. Inicializar tripodpy ───────────────────────────────────────────────
    pipeline.initialize_simulation()

    # ── 3. Química del disco ─────────────────────────────────────────────────
    pipeline.add_volatile_components()

    # ── 4. Física de snowlines ───────────────────────────────────────────────
    pipeline.setup_physics()
    pipeline.setup_star_evolution()
    pipeline.add_snowline_fields()

    # ── 5. Gap planetario (si aplica) ────────────────────────────────────────
    gap_model = cfg.get("gap_model", "none")
    if gap_model != "none":
        pipeline.gap_alpha_ref   = cfg["alpha_ref"]
        pipeline.gap_M_planet    = cfg["M_planet"]
        pipeline.gap_a_planet_au = cfg["a_planet_au"]

        if gap_model == "kanagawa":
            pipeline.setup_gap_kanagawa(imprint=cfg.get("imprint", False))
        elif gap_model == "duffell":
            pipeline.setup_gap_duffell(imprint=cfg.get("imprint", False))
        else:
            print(f"  [!] gap_model='{gap_model}' desconocido. Sin gap.")

    # ── 6. Fields de composición → HDF5 ──────────────────────────────────────
    pipeline.add_ice_sigma_fields()

    # ── 7. Sincronización y run ───────────────────────────────────────────────
    pipeline.sim.update()
    pipeline.run_integration(
        t_end_years   = cfg["t_end_years"],
        num_snapshots = cfg["num_snapshots"],
    )

    elapsed = time.time() - t0
    print(f"\n  ✓ Run '{label}' completado en {elapsed / 60:.1f} min")
    return datadir, elapsed


# ══════════════════════════════════════════════════════════════════════════════
# Función: generar diagnósticos de UN run
# ══════════════════════════════════════════════════════════════════════════════

def run_diagnostics_single(cfg: dict, show: bool = False):
    """
    Genera todas las figuras de diagnóstico para un run individual.

    Parameters
    ----------
    cfg : dict
        Configuración del run.
    show : bool
        Si True, llama plt.show() al final.
    """
    label   = cfg["label"]
    datadir = cfg["datadir"]
    savedir = os.path.join(datadir, "figures")

    print(f"\n  → Diagnósticos: {label}  →  {savedir}")

    try:
        d = SnowlineDiagnostics(datadir, savedir=savedir, r_trim=0.93, t_min_yr=100.0)
    except Exception as e:
        print(f"  [!] Error cargando {datadir}: {e}")
        return

    # Figuras estándar del pipeline
    try: d.plot_hovmoller(quantity="eps")
    except Exception as e: print(f"     hovmoller_eps: {e}")

    try: d.plot_hovmoller(quantity="Sigma_dust")
    except Exception as e: print(f"     hovmoller_Sigma_dust: {e}")

    try: d.plot_hovmoller_rt(quantities=["Sigma_gas", "T", "a_max"], t_unit="kyr")
    except Exception as e: print(f"     hovmoller_rt: {e}")

    try: d.plot_size_distribution(it=-1)
    except Exception as e: print(f"     size_distribution: {e}")

    try: d.plot_pebble_flux()
    except Exception as e: print(f"     pebble_flux: {e}")

    try: d.plot_profiles(it=-1)
    except Exception as e: print(f"     profiles: {e}")

    try: d.plot_hovmoller_comp(quantity='dust')
    except Exception as e: print(f"     hovmoller_comp_dust: {e}")

    try: d.plot_hovmoller_comp(quantity='gas')
    except Exception as e: print(f"     hovmoller_comp_gas: {e}")

    if show:
        plt.show()
    plt.close("all")
    print(f"Figuras guardadas en: {savedir}")


# ══════════════════════════════════════════════════════════════════════════════
# Función: figura comparativa entre todos los runs
# ══════════════════════════════════════════════════════════════════════════════

def plot_comparison(runs: list, savedir: str = FIGURES_DIR, show: bool = False):
    """
    Genera figuras de comparación entre todos los runs activos.

    Panel 1: Σ_dust(r) en el snapshot final — superpuestos
    Panel 2: Σ_gas(r) en el snapshot final — superpuestos
    Panel 3: a_max(r) en el snapshot final — superpuestos
    Panel 4: rsnow_H2O(t) — evolución temporal del snowline — superpuestos
    """
    print(f"\n{'─'*60}")
    print("  Generando figura de comparación multi-run...")
    os.makedirs(savedir, exist_ok=True)

    # Colores para cada run
    PALETTE = ["#90CAF9", "#EF9A9A", "#A5D6A7", "#FFCC80", "#CE93D8"]
    LINESTYLES = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]

    fig = plt.figure(figsize=(14, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)
    ax_dust  = fig.add_subplot(gs[0, 0])
    ax_gas   = fig.add_subplot(gs[0, 1])
    ax_amax  = fig.add_subplot(gs[1, 0])
    ax_rsnow = fig.add_subplot(gs[1, 1])

    for idx, cfg in enumerate(runs):
        label   = cfg["label"]
        datadir = cfg["datadir"]
        color   = PALETTE[idx % len(PALETTE)]
        ls      = LINESTYLES[idx % len(LINESTYLES)]

        # Verificar que el datadir exista
        if not os.path.isdir(datadir):
            print(f"  [!] No encontrado: {datadir}  (¿corriste el run?)")
            continue

        try:
            d = SnowlineDiagnostics(datadir, savedir=None, r_trim=0.93, t_min_yr=100.0)
        except Exception as e:
            print(f"  [!] Error leyendo {label}: {e}")
            continue

        r_au = d.r_au
        it   = d.Nt - 1   # último snapshot

        # ── Panel 1: Σ_dust ───────────────────────────────────────────────────
        try:
            Sig_d = d.data.dust.Sigma[np.where(d.t_mask)[0][it], d.r_mask, :].sum(-1)
            ax_dust.loglog(r_au, Sig_d, color=color, ls=ls, lw=1.8, label=label)
        except Exception: pass

        # ── Panel 2: Σ_gas ────────────────────────────────────────────────────
        try:
            Sig_g = d.data.gas.Sigma[np.where(d.t_mask)[0][it], d.r_mask]
            ax_gas.loglog(r_au, Sig_g, color=color, ls=ls, lw=1.8, label=label)
        except Exception: pass

        # ── Panel 3: a_max ────────────────────────────────────────────────────
        try:
            a_max = d.data.dust.s.max[np.where(d.t_mask)[0][it], d.r_mask]
            ax_amax.loglog(r_au, a_max, color=color, ls=ls, lw=1.8, label=label)
        except Exception: pass

        # ── Panel 4: rsnow_H2O(t) ─────────────────────────────────────────────
        try:
            r_snow = d._get_rsnow_series("H2O")
            if r_snow is not None:
                ax_rsnow.plot(d.t_yr, r_snow, color=color, ls=ls, lw=1.8, label=label)
        except Exception: pass

    # Decorar paneles
    t_last = cfg.get("t_end_years", 1e6)   # último run como referencia

    for ax, ylabel, title in [
        (ax_dust,  r"$\Sigma_{dust}$ [g cm$^{-2}$]", f"Σ_dust  |  t = {t_last:.0e} yr"),
        (ax_gas,   r"$\Sigma_{gas}$ [g cm$^{-2}$]",  f"Σ_gas   |  t = {t_last:.0e} yr"),
        (ax_amax,  r"$a_{max}$ [cm]",                 f"a_max   |  t = {t_last:.0e} yr"),
    ]:
        ax.set_xlabel("Radio [AU]")
        ax.set_ylabel(ylabel)
        ax.set_title(title, pad=6)
        ax.legend(fontsize=7, framealpha=0.5)

    ax_rsnow.set_xscale("log")
    ax_rsnow.set_xlabel("Tiempo [yr]")
    ax_rsnow.set_ylabel(r"$r_{snow,H_2O}$ [AU]")
    ax_rsnow.set_title(r"Evolución del snowline de H$_2$O", pad=6)
    ax_rsnow.legend(fontsize=7, framealpha=0.5)

    fig.suptitle("Comparación multi-run  —  PPOLs Pipeline",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()

    out_path = os.path.join(savedir, "comparison_all_runs.pdf")
    fig.savefig(out_path, bbox_inches="tight")
    print(f"  ✓ Figura comparativa guardada: {out_path}")

    if show:
        plt.show()
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# Punto de entrada
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="PPOLs multi-run: simular + diagnósticos + comparación"
    )
    parser.add_argument(
        "--diag", action="store_true",
        help="Solo genera diagnósticos (no corre simulaciones de nuevo)"
    )
    parser.add_argument(
        "--compare", action="store_true",
        help="Solo genera la figura comparativa"
    )
    parser.add_argument(
        "--runs", nargs="+", metavar="LABEL",
        help="Labels específicos a procesar (default: todos)"
    )
    parser.add_argument(
        "--show", action="store_true",
        help="Llamar plt.show() al final de cada bloque de figuras"
    )
    args = parser.parse_args()

    # Filtrar runs por label si se especificaron
    selected_runs = RUNS
    if args.runs:
        selected_runs = [r for r in RUNS if r["label"] in args.runs]
        if not selected_runs:
            print(f"[!] Ningún run coincide con: {args.runs}")
            print(f"    Labels disponibles: {[r['label'] for r in RUNS]}")
            sys.exit(1)

    os.makedirs(BASE_DIR, exist_ok=True)

    # ── Fase 1: Simulaciones ─────────────────────────────────────────────────
    if not args.diag and not args.compare:
        print("\n" + "═"*60)
        print(f"  FASE 1: Corriendo {len(selected_runs)} run(s)")
        print("═"*60)

        total_t0 = time.time()
        for cfg in selected_runs:
            try:
                build_and_run(cfg)
            except Exception as e:
                print(f"\n  [ERROR] Run '{cfg['label']}': {e}")
                import traceback; traceback.print_exc()

        total_min = (time.time() - total_t0) / 60
        print(f"\n  ✓ Todas las simulaciones completadas en {total_min:.1f} min")

    # ── Fase 2: Diagnósticos individuales ────────────────────────────────────
    if not args.compare:
        print("\n" + "═"*60)
        print("  FASE 2: Generando diagnósticos individuales")
        print("═"*60)

        for cfg in selected_runs:
            run_diagnostics_single(cfg, show=args.show)

    # ── Fase 3: Figura comparativa ───────────────────────────────────────────
    print("\n" + "═"*60)
    print("  FASE 3: Figura comparativa multi-run")
    print("═"*60)

    plot_comparison(selected_runs, savedir=FIGURES_DIR, show=args.show)

    print("\n" + "═"*60)
    print("  run_diagnostics.py completado.")
    print(f"  Figuras en: {FIGURES_DIR}")
    print("═"*60)


if __name__ == "__main__":
    main()

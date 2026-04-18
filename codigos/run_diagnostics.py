"""
run_diagnostics.py — Script de post-procesamiento multi-run (gaps & substructures)
====================================================================================

Objetivo científico
-------------------
Estudiar cómo los gaps y sub-estructuras del disco protoplanetario afectan
la composición de protoplanetas acretados.

Todas las simulaciones usan:
  - Modelo de gap: Duffell (2020) — preferido por su perfil empírico más realista
  - Especie activa: H2O (snowline principal)
  - Integración: 10 Myr, 60 snapshots log-uniformes
  - Estrella: 1 M☉
  - Grilla: 200 celdas, 1–300 AU

Grupos de configuraciones
--------------------------
1. MISMO TAMAÑO, DISTINTA POSICIÓN
   Un gap de ~Júpiter a 2, 5, 10 y 20 AU

2. MISMA POSICIÓN (5 AU), DISTINTO TAMAÑO
   Gaps a 5 AU con 0.1, 0.5, 1.0 y 2.0 M_Júpiter

3. MÚLTIPLES GAPS EN POSICIONES DISTINTAS
   3a: Júpiter@5AU,  (½MJ)@7AU,  Júpiter@10AU
   3b: Saturno@5AU,  (½MJ)@7AU,  Júpiter@10AU

4. GAPS SINUSOIDALES (múltiples gaps equidistantes)
   Amplitud suave (A=1), media (A=5) y fuerte (A=10)

Uso
---
    python run_diagnostics.py            → simula + diagnósticos + comparación
    python run_diagnostics.py --diag     → solo diagnósticos
    python run_diagnostics.py --compare  → solo figura comparativa
    python run_diagnostics.py --runs baseline duffell_5au_1MJ
    python run_diagnostics.py --show     → plt.show() interactivo
"""

import pipeline_snowlines
import os
import sys
import time
import argparse

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import dustpy.constants as c

from pipeline_snowlines import WaterworldPipeline
from plot_diagnostics    import SnowlineDiagnostics

# ══════════════════════════════════════════════════════════════════════════════
# Directorios de salida
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR    = "data_gaps_pipeline/t_1e5"
FIGURES_DIR = "figs_gaps_pipeline/t_1e5"

# ══════════════════════════════════════════════════════════════════════════════
# Parámetros globales compartidos por todos los runs
# ══════════════════════════════════════════════════════════════════════════════
_DEFAULTS = dict(
    active_species = ["H2O"],
    t_end_years    = 1e5,     # 5 Myr
    num_snapshots  = 50,
    M_star_Msun    = 1.0,
    Nr             = 100,      # Nr=200 colapsa el solver implícito (N≈28k → muy lento)
    rmax_au        = 200.0,    # 300→200 AU: el polvo relevante está <100 AU
    alpha_ref      = 1e-3,
)

def _run(label, gap_model, **kwargs):
    """Helper para construir un dict de configuración completo."""
    cfg = dict(_DEFAULTS)
    cfg.update(kwargs)
    cfg["label"]     = label
    cfg["datadir"]   = os.path.join(BASE_DIR, label)
    cfg["gap_model"] = gap_model
    cfg.setdefault("M_planet",    0.0)
    cfg.setdefault("a_planet_au", 0.0)
    cfg.setdefault("imprint",     False)
    # Sinusoidal params (solo para gap_model="sinusoidal")
    cfg.setdefault("amplitude",   5.0)
    cfg.setdefault("n_bumps",     5)
    cfg.setdefault("r_inner_au",  1.0)
    cfg.setdefault("r_outer_au",  100.0)
    return cfg


# ══════════════════════════════════════════════════════════════════════════════
# Configuraciones de runs
# ══════════════════════════════════════════════════════════════════════════════

RUNS = [

    # ── Baseline ─────────────────────────────────────────────────────────────
    _run("baseline", "none"),

    # ── Grupo 1: Mismo tamaño (1 M_Júpiter), distintas posiciones ────────────
    _run("duffell_1MJ_2au",  "duffell", M_planet=317.8*c.M_earth, a_planet_au=2.0),
    _run("duffell_1MJ_5au",  "duffell", M_planet=317.8*c.M_earth, a_planet_au=5.0),
    _run("duffell_1MJ_10au", "duffell", M_planet=317.8*c.M_earth, a_planet_au=10.0),
    _run("duffell_1MJ_20au", "duffell", M_planet=317.8*c.M_earth, a_planet_au=20.0),

    # ── Grupo 2: Misma posición (5 AU), distinto tamaño ──────────────────────
    _run("duffell_01MJ_5au", "duffell", M_planet=0.1*317.8*c.M_earth, a_planet_au=5.0),
    _run("duffell_05MJ_5au", "duffell", M_planet=0.5*317.8*c.M_earth, a_planet_au=5.0),
    # duffell_1MJ_5au ya definido arriba
    _run("duffell_2MJ_5au",  "duffell", M_planet=2.0*317.8*c.M_earth, a_planet_au=5.0),

    # ── Grupo 3a: Múltiples gaps — 3 planetas, distintas posiciones ──────────
    # Júpiter@5 + ½MJ@7 + Júpiter@10 AU
    _run("duffell_multi_3a_5_7_10",  "duffell_multi",
         planets=[
             {"M_planet": 317.8*c.M_earth,       "a_planet_au": 5.0},
             {"M_planet": 0.5*317.8*c.M_earth,   "a_planet_au": 7.0},
             {"M_planet": 317.8*c.M_earth,        "a_planet_au": 10.0},
         ]),

    # ── Grupo 3b: Saturno@5 + ½MJ@7 + Júpiter@10 AU ─────────────────────────
    _run("duffell_multi_3b_sat_half_jup", "duffell_multi",
         planets=[
             {"M_planet": 95.159*c.M_earth,       "a_planet_au": 5.0},   # Saturno
             {"M_planet": 0.5*317.8*c.M_earth,    "a_planet_au": 7.0},   # ½MJ
             {"M_planet": 317.8*c.M_earth,         "a_planet_au": 10.0},  # Júpiter
         ]),

    # ── Grupo 4: Gaps sinusoidales ────────────────────────────────────────────
    _run("sinusoidal_A1_suave",  "sinusoidal", amplitude=1.0,  n_bumps=5,
         r_inner_au=5.0, r_outer_au=100.0),
    _run("sinusoidal_A5_medio",  "sinusoidal", amplitude=5.0,  n_bumps=5,
         r_inner_au=5.0, r_outer_au=100.0),
    _run("sinusoidal_A10_fuerte","sinusoidal", amplitude=10.0, n_bumps=5,
         r_inner_au=5.0, r_outer_au=100.0),
]


# ══════════════════════════════════════════════════════════════════════════════
# Construcción y ejecución de un run
# ══════════════════════════════════════════════════════════════════════════════

def build_and_run(cfg: dict):
    label   = cfg["label"]
    datadir = cfg["datadir"]

    sep = "═" * 60
    print(f"\n{sep}\n  RUN: {label}\n  datadir: {datadir}\n{sep}")

    t0       = time.time()
    pipeline = WaterworldPipeline(datadir=datadir)

    # 1. Grilla y estrella
    pipeline.active_species = cfg["active_species"]
    pipeline.setup_grid(rmin=1.0*c.au, rmax=cfg["rmax_au"]*c.au, Nr=cfg["Nr"])
    pipeline.setup_star(M_star_Msun=cfg["M_star_Msun"])

    # 2. Inicializar tripodpy
    pipeline.initialize_simulation()

    # 3. Química
    pipeline.add_volatile_components()

    # 4. Snowlines
    pipeline.setup_physics()
    pipeline.setup_star_evolution()
    pipeline.add_snowline_fields()

    # 5. Gap / sub-estructura
    gap_model = cfg.get("gap_model", "none")
    pipeline.gap_alpha_ref = cfg["alpha_ref"]

    if gap_model == "duffell":
        pipeline.gap_M_planet    = cfg["M_planet"]
        pipeline.gap_a_planet_au = cfg["a_planet_au"]
        pipeline.setup_gap_duffell(imprint=cfg.get("imprint", False))

    elif gap_model == "duffell_multi":
        # Múltiples gaps: aplicamos Duffell secuencialmente sobreescribiendo α
        # Nota: cada llamada sobreescribe el updater anterior, por lo que
        # instalamos un updater compuesto que aplica TODOS los planetas
        _planets  = cfg["planets"]
        _a0_val   = cfg["alpha_ref"]
        alpha0    = pipeline.sim.gas.alpha.copy()

        from scipy.interpolate import interp1d
        from pipeline_methods.pressure_bumps import duffell2020

        _pl_list  = [(p["M_planet"], p["a_planet_au"] * c.au) for p in _planets]
        _a0_cl    = _a0_val
        _alpha0   = alpha0

        print(f"  → {len(_pl_list)} gaps Duffell (multi):")
        for Mp, a in _pl_list:
            print(f"     M={Mp/c.M_earth:.1f} M⊕  a={a/c.au:.1f} AU  "
                  f"q={Mp/pipeline.sim.star.M:.2e}")

        def _multi_duffell(sim):
            alpha_out = _alpha0.copy()
            for Mp, a_pl in _pl_list:
                q      = Mp / float(sim.star.M)
                f_h    = interp1d(sim.grid.r, sim.gas.Hp / sim.grid.r,
                                  bounds_error=False, fill_value="extrapolate")
                h_p    = float(f_h(a_pl))
                f_alp  = interp1d(sim.grid.r, _alpha0,
                                  bounds_error=False, fill_value="extrapolate")
                alp_p  = float(f_alp(a_pl))
                f_gap  = duffell2020(sim.grid.r, a_pl, q, h_p, alp_p)
                f_safe = np.maximum(f_gap, 1e-10)
                alpha_out /= f_safe
            return alpha_out

        pipeline.sim.gas.alpha.updater.updater = _multi_duffell
        pipeline.sim.gas.alpha.update()
        pipeline.sim.update()

    elif gap_model == "sinusoidal":
        pipeline.setup_alpha_sinusoidal(
            alpha_ref  = cfg["alpha_ref"],
            amplitude  = cfg["amplitude"],
            n_bumps    = cfg["n_bumps"],
            r_inner_au = cfg["r_inner_au"],
            r_outer_au = cfg["r_outer_au"],
            imprint    = cfg.get("imprint", False),
        )

    elif gap_model != "none":
        print(f"  [!] gap_model='{gap_model}' desconocido — corriendo sin gap.")

    # 6. Fields HDF5
    pipeline.add_ice_sigma_fields()

    # 7. Run
    pipeline.sim.update()
    pipeline.run_integration(
        t_start_years = 1,
        t_end_years   = cfg["t_end_years"],
        num_snapshots = cfg["num_snapshots"],
    )

    elapsed = time.time() - t0
    print(f"\n  ✓ '{label}' completado en {elapsed/60:.1f} min")
    return datadir, elapsed


# ══════════════════════════════════════════════════════════════════════════════
# Diagnósticos individuales
# ══════════════════════════════════════════════════════════════════════════════

def run_diagnostics_single(cfg: dict, show: bool = False):
    label   = cfg["label"]
    datadir = cfg["datadir"]
    savedir = os.path.join(datadir, "figures")

    print(f"\n  → Diagnósticos: {label}  →  {savedir}")
    try:
        d = SnowlineDiagnostics(datadir, savedir=savedir, r_trim=0.93, t_min_yr=100.0)
    except Exception as e:
        print(f"  [!] Error cargando {datadir}: {e}")
        return

    for method, kwargs in [
        ("plot_hovmoller",       {"quantity": "eps"}),
        ("plot_hovmoller",       {"quantity": "Sigma_dust"}),
        ("plot_hovmoller_rt",    {"quantities": ["Sigma_gas", "T", "a_max"], "t_unit": "Myr"}),
        ("plot_size_distribution", {"it": -1}),
        ("plot_pebble_flux",     {}),
        ("plot_profiles",        {"it": -1}),
        ("plot_hovmoller_comp",  {"quantity": "dust"}),
        ("plot_hovmoller_comp",  {"quantity": "gas"}),
    ]:
        try:
            getattr(d, method)(**kwargs)
        except Exception as e:
            print(f"     {method}: {e}")

    if show:
        plt.show()
    plt.close("all")
    print(f"Figuras guardadas en: {savedir}")


# ══════════════════════════════════════════════════════════════════════════════
# Figura comparativa global
# ══════════════════════════════════════════════════════════════════════════════

def plot_comparison(runs: list, savedir: str = FIGURES_DIR, show: bool = False):
    print(f"\n{'─'*60}\n  Generando figura de comparación multi-run...")
    os.makedirs(savedir, exist_ok=True)

    PALETTE  = plt.cm.tab20(np.linspace(0, 1, len(runs)))
    LINESTYLES = ["-", "--", "-.", ":", (0, (3,1,1,1)), (0, (5,1))]

    fig = plt.figure(figsize=(16, 10))
    gs  = gridspec.GridSpec(2, 2, figure=fig, hspace=0.38, wspace=0.32)
    ax_dust  = fig.add_subplot(gs[0, 0])
    ax_gas   = fig.add_subplot(gs[0, 1])
    ax_amax  = fig.add_subplot(gs[1, 0])
    ax_rsnow = fig.add_subplot(gs[1, 1])

    for idx, cfg in enumerate(runs):
        label   = cfg["label"]
        datadir = cfg["datadir"]
        color   = PALETTE[idx]
        ls      = LINESTYLES[idx % len(LINESTYLES)]

        if not os.path.isdir(datadir):
            print(f"  [!] No encontrado: {datadir}")
            continue

        try:
            d = SnowlineDiagnostics(datadir, savedir=None, r_trim=0.93, t_min_yr=100.0)
        except Exception as e:
            print(f"  [!] Error leyendo {label}: {e}")
            continue

        r_au = d.r_au
        it   = d.Nt - 1
        it_f = np.where(d.t_mask)[0][it]

        for ax, field, label_str in [
            (ax_dust, lambda: d.data.dust.Sigma[it_f, d.r_mask, :].sum(-1), "Σ_dust"),
            (ax_gas,  lambda: d.data.gas.Sigma[it_f, d.r_mask],             "Σ_gas"),
            (ax_amax, lambda: d.data.dust.s.max[it_f, d.r_mask],            "a_max"),
        ]:
            try:
                ax.loglog(r_au, field(), color=color, ls=ls, lw=1.6, label=label)
            except Exception:
                pass

        try:
            r_snow = d._get_rsnow_series("H2O")
            if r_snow is not None:
                ax_rsnow.plot(d.t_yr, r_snow, color=color, ls=ls, lw=1.6, label=label)
        except Exception:
            pass

    t_last = _DEFAULTS["t_end_years"]
    for ax, ylabel, title in [
        (ax_dust,  r"$\Sigma_{dust}$ [g cm$^{-2}$]", f"Σ_dust  |  t = {t_last:.0e} yr"),
        (ax_gas,   r"$\Sigma_{gas}$ [g cm$^{-2}$]",  f"Σ_gas   |  t = {t_last:.0e} yr"),
        (ax_amax,  r"$a_{max}$ [cm]",                 f"a_max   |  t = {t_last:.0e} yr"),
    ]:
        ax.set_xlabel("Radio [AU]"); ax.set_ylabel(ylabel); ax.set_title(title, pad=6)
        ax.legend(fontsize=6, framealpha=0.5, ncol=2)

    ax_rsnow.set_xscale("log")
    ax_rsnow.set_xlabel("Tiempo [yr]")
    ax_rsnow.set_ylabel(r"$r_\mathrm{snow,H_2O}$ [AU]")
    ax_rsnow.set_title(r"Snowline H$_2$O", pad=6)
    ax_rsnow.legend(fontsize=6, framealpha=0.5, ncol=2)

    fig.suptitle("Comparación multi-run — Efecto de gaps en composición de protoplanetas",
                 fontsize=13, fontweight="bold", y=1.01)
    fig.tight_layout()
    out = os.path.join(savedir, "comparison_all_runs.pdf")
    fig.savefig(out, bbox_inches="tight")
    print(f"  ✓ Guardado: {out}")
    if show: plt.show()
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════════
# CLI
# ══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="PPOLs multi-run: gaps & sub-estructuras → composición de protoplanetas"
    )
    parser.add_argument("--diag",    action="store_true", help="Solo diagnósticos")
    parser.add_argument("--compare", action="store_true", help="Solo figura comparativa")
    parser.add_argument("--runs", nargs="+", metavar="LABEL",
                        help="Labels específicos (default: todos)")
    parser.add_argument("--show",    action="store_true", help="plt.show() interactivo")
    args = parser.parse_args()

    selected = RUNS
    if args.runs:
        selected = [r for r in RUNS if r["label"] in args.runs]
        if not selected:
            print(f"[!] Labels disponibles: {[r['label'] for r in RUNS]}")
            sys.exit(1)

    os.makedirs(BASE_DIR, exist_ok=True)

    if not args.diag and not args.compare:
        print(f"\n{'═'*60}\n  FASE 1: {len(selected)} simulación(es) — {_DEFAULTS['t_end_years']:.0e} yr cada una\n{'═'*60}")
        t0 = time.time()
        for cfg in selected:
            try: build_and_run(cfg)
            except Exception as e:
                import traceback
                print(f"\n  [ERROR] '{cfg['label']}': {e}")
                traceback.print_exc()
        print(f"\n  ✓ Simulaciones completadas en {(time.time()-t0)/60:.1f} min")

    if not args.compare:
        print(f"\n{'═'*60}\n  FASE 2: Diagnósticos individuales\n{'═'*60}")
        for cfg in selected:
            run_diagnostics_single(cfg, show=args.show)

    print(f"\n{'═'*60}\n  FASE 3: Comparación global\n{'═'*60}")
    plot_comparison(selected, savedir=FIGURES_DIR, show=args.show)

    print(f"\n{'═'*60}\n  Completado. Figuras → {FIGURES_DIR}/\n{'═'*60}")


if __name__ == "__main__":
    main()

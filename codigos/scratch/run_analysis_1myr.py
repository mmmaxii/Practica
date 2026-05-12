"""
run_analysis_1myr.py
====================
Análisis científico de las runs de data_1myr/ con PA3.

Figuras generadas:
  1. heatmap_fh2o.pdf      — gap_pos × gap_mass → f_H2O final  (solo single_*)
  2. bar_fh2o.pdf          — todas las runs ordenadas por f_H2O final
  3. mass_grouped.pdf      — M(t) agrupado por categoría (3 paneles)
  4. pebble_flux_snow.pdf  — flujo de pebbles en la snowline

Uso:
    py run_analysis_1myr.py --r_planet 2.0 --savedir figs_1myr
"""

import os, sys, re, glob, argparse

# Forzar UTF-8 en consola Windows
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import dustpy.constants as c

_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_HERE, "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3

# ── Constantes ────────────────────────────────────────────────────────────────
AU     = c.au
YR     = c.year
M_EARTH = c.M_earth

# ── Estética global oscura ────────────────────────────────────────────────────
BG  = "#0d1117"
FG  = "#e6edf3"
GRID_C = "#21262d"

plt.rcParams.update({
    "figure.facecolor":    BG,
    "axes.facecolor":      "#161b22",
    "axes.edgecolor":      GRID_C,
    "axes.labelcolor":     FG,
    "axes.titlecolor":     FG,
    "xtick.color":         FG,
    "ytick.color":         FG,
    "text.color":          FG,
    "grid.color":          GRID_C,
    "grid.linewidth":      0.6,
    "axes.grid":           True,
    "legend.facecolor":    "#161b22",
    "legend.edgecolor":    GRID_C,
    "font.size":           10,
    "axes.titlesize":      13,
    "axes.labelsize":      11,
    "savefig.facecolor":   BG,
    "savefig.dpi":         160,
})

# ── Mapeo de masa del gap ─────────────────────────────────────────────────────
MASS_LABELS = {
    "sup_earth": "Super-Earth",
    "nep":       "Neptune",
    "sat":       "Saturn",
    "jup":       "Jupiter",
    "sup_jup":   "Super-Jupiter",
}
MASS_ORDER = ["sup_earth", "nep", "sat", "jup", "sup_jup"]

# Colores de composición
C_SIL = "#ef9a9a"
C_CO2 = "#80cbc4"
C_H2O = "#4fc3f7"

# ── Utilidades ────────────────────────────────────────────────────────────────

def collect_runs(dataroot):
    runs = []
    for entry in sorted(os.listdir(dataroot)):
        path = os.path.join(dataroot, entry)
        if os.path.isdir(path) and glob.glob(os.path.join(path, "data*.hdf5")):
            runs.append((entry, path))
    return runs


def parse_single(name):
    """'single_jup_3.0au' → {'mass':'jup', 'pos':3.0} o None"""
    m = re.match(r"single_([a-z_]+)_([\d.]+)au", name)
    if not m:
        return None
    mass_raw = m.group(1)
    # normalizar sup_earth y sup_jup
    if "sup_earth" in mass_raw:
        mass = "sup_earth"
    elif "sup_jup" in mass_raw:
        mass = "sup_jup"
    else:
        mass = mass_raw
    return {"mass": mass, "pos": float(m.group(2))}


def group_runs(runs):
    single, multi, sinus = [], [], []
    for name, path in runs:
        if name.startswith("single_"):
            meta = parse_single(name)
            single.append({"name": name, "path": path, "meta": meta})
        elif name.startswith("multi_"):
            multi.append({"name": name, "path": path, "meta": None})
        elif name.startswith("sinusoidal_"):
            sinus.append({"name": name, "path": path, "meta": None})
    return single, multi, sinus


def run_pa3(path, r_planet_au, t_min_yr=100.0):
    pa3  = PebbleAccretionModule3.from_datadir(path, t_min_yr=t_min_yr)
    hist = pa3.run_growth([r_planet_au])[r_planet_au]
    return pa3, hist


def final_fh2o(hist):
    """Fracción final de H2O en %"""
    if len(hist) == 0:
        return 0.0
    M = hist[-1, 1]
    H = hist[-1, 2]
    return 100.0 * H / (M + 1e-30)


def final_mass_earth(hist):
    if len(hist) == 0:
        return 0.0
    return hist[-1, 1] / M_EARTH


def pretty(name):
    return name.replace("_", " ")


# ═══════════════════════════════════════════════════════════════════════════════
# Figura 1: Heatmap gap_pos × gap_mass → f_H2O   (solo single_*)
# ═══════════════════════════════════════════════════════════════════════════════

def plot_heatmap(single_data, savedir, r_planet_au):
    """single_data: list de dicts con keys name, path, meta, hist"""
    positions = sorted({d["meta"]["pos"] for d in single_data if d["meta"]})
    masses    = MASS_ORDER

    # Construir matriz
    Z = np.full((len(masses), len(positions)), np.nan)
    M_mat = np.full_like(Z, np.nan)
    for d in single_data:
        if not d["meta"] or "hist" not in d:
            continue
        mi = masses.index(d["meta"]["mass"]) if d["meta"]["mass"] in masses else None
        pi = positions.index(d["meta"]["pos"]) if d["meta"]["pos"] in positions else None
        if mi is None or pi is None:
            continue
        Z[mi, pi]     = final_fh2o(d["hist"])
        M_mat[mi, pi] = final_mass_earth(d["hist"])

    fig, ax = plt.subplots(figsize=(10, 5))

    cmap = cm.get_cmap("Blues")
    cmap.set_bad("#1a1f2e")
    im = ax.imshow(Z, cmap=cmap, aspect="auto", origin="lower",
                   vmin=0, vmax=max(15, np.nanmax(Z)))

    # Etiquetas en cada celda
    for mi in range(len(masses)):
        for pi in range(len(positions)):
            val = Z[mi, pi]
            mval = M_mat[mi, pi]
            if not np.isnan(val):
                ww = val > 10
                color = "#0d1117" if val > 6 else FG
                ax.text(pi, mi, f"{val:.1f}%\n{mval*1000:.2f}m⊕",
                        ha="center", va="center", fontsize=7.5,
                        color=color, fontweight="bold" if ww else "normal")
                if ww:
                    ax.add_patch(plt.Rectangle(
                        (pi - 0.5, mi - 0.5), 1, 1,
                        fill=False, edgecolor="#ffd700", lw=2.5, zorder=5))

    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("$f_{\\mathrm{H_2O}}$ final [%]", color=FG, fontsize=10)
    cb.ax.tick_params(colors=FG)
    cb.outline.set_edgecolor(GRID_C)

    ax.set_xticks(range(len(positions)))
    ax.set_xticklabels([f"{p:.0f} AU" for p in positions])
    ax.set_yticks(range(len(masses)))
    ax.set_yticklabels([MASS_LABELS.get(m, m) for m in masses])
    ax.set_xlabel("Posición del gap [AU]")
    ax.set_ylabel("Masa del planeta del gap")
    ax.set_title(f"Fracción final de H₂O  —  embrión @ {r_planet_au:.1f} AU\n"
                 f"(borde dorado = waterworld > 10%)", pad=10)
    ax.grid(False)

    fig.tight_layout()
    out = os.path.join(savedir, "heatmap_fh2o.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"  [OK] {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# Figura 2: Bar chart ordenado por f_H2O — todas las runs
# ═══════════════════════════════════════════════════════════════════════════════

def plot_bar_fh2o(all_data, savedir, r_planet_au):
    """all_data: list de dicts con keys name, hist, group"""
    entries = [(d["name"], final_fh2o(d["hist"]), final_mass_earth(d["hist"]), d["group"])
               for d in all_data if "hist" in d and len(d["hist"]) > 0]
    entries.sort(key=lambda x: x[1], reverse=True)

    if not entries:
        return

    names, fh2o_vals, mass_vals, groups = zip(*entries)
    n = len(names)
    y = np.arange(n)

    group_color = {"single": "#4fc3f7", "multi": "#ff8a65", "sinusoidal": "#a5d6a7"}

    fig, ax = plt.subplots(figsize=(10, max(5, n * 0.28 + 2)))

    bars = ax.barh(y, fh2o_vals, height=0.7, zorder=3,
                   color=[group_color.get(g, "#aaa") for g in groups])

    # Línea de waterworld
    ax.axvline(10, color="#ffd700", lw=1.5, ls="--", alpha=0.8, zorder=4,
               label="Límite waterworld (10%)")

    # Etiqueta de masa final
    for i, (fh2o, mass) in enumerate(zip(fh2o_vals, mass_vals)):
        ax.text(fh2o + 0.1, i, f"{mass*1000:.2f} m⊕",
                va="center", ha="left", fontsize=7, color=FG, alpha=0.8)

    ax.set_yticks(y)
    ax.set_yticklabels([pretty(n) for n in names], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlabel("$f_{\\mathrm{H_2O}}$ final [%]")
    ax.set_title(f"Fracción de agua final por run  —  embrión @ {r_planet_au:.1f} AU", pad=10)

    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=v, label=k) for k, v in group_color.items()]
    legend_els.append(plt.Line2D([0], [0], color="#ffd700", ls="--", label="WW threshold 10%"))
    ax.legend(handles=legend_els, loc="lower right", fontsize=8, framealpha=0.4)

    fig.tight_layout()
    out = os.path.join(savedir, "bar_fh2o.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"  [OK] {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# Figura 3: Evolución temporal agrupada — 3 paneles
# ═══════════════════════════════════════════════════════════════════════════════

def plot_mass_grouped(groups_data, savedir, r_planet_au):
    """groups_data: dict con keys 'single', 'multi', 'sinusoidal', cada uno lista de dicts"""
    titles = {"single": "Single planet gap", "multi": "Multi-planet gaps",
              "sinusoidal": "Sinusoidal α profile"}
    group_keys = [k for k in ["single", "multi", "sinusoidal"] if groups_data.get(k)]
    n_panels = len(group_keys)

    fig, axes = plt.subplots(1, n_panels, figsize=(5.5 * n_panels, 5.5),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, gkey in zip(axes, group_keys):
        data = groups_data[gkey]
        # Colormap por posición o índice
        cmap = cm.get_cmap("plasma", max(len(data), 1))

        for idx, d in enumerate(data):
            hist = d.get("hist")
            if hist is None or len(hist) == 0:
                continue
            t_yr = hist[:, 0] / YR
            M_Me = hist[:, 1] / M_EARTH
            fh2o = 100 * hist[:, 2] / (hist[:, 1] + 1e-30)

            color = cmap(idx / max(len(data) - 1, 1))
            lw = 1.5
            if gkey == "single" and d.get("meta"):
                # Para single: color por posición del gap
                positions = sorted({dd["meta"]["pos"] for dd in data if dd.get("meta")})
                pos = d["meta"]["pos"]
                cidx = positions.index(pos) / max(len(positions) - 1, 1)
                color = cm.get_cmap("cool")(cidx)
                # linestyle por masa
                mass_ls = {"sup_earth": ":", "nep": "-.", "sat": "--",
                           "jup": "-", "sup_jup": (0, (5, 1))}
                ls = mass_ls.get(d["meta"]["mass"], "-")
                ax.plot(t_yr, M_Me, color=color, lw=lw, ls=ls, alpha=0.85)
            else:
                ax.plot(t_yr, M_Me, color=color, lw=lw, alpha=0.85,
                        label=pretty(d["name"]))

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Tiempo [años]")
        ax.set_title(titles[gkey])
        if ax == axes[0]:
            ax.set_ylabel(r"$M_{\rm core}$ [$M_\oplus$]")

        if gkey != "single":
            ax.legend(fontsize=7, framealpha=0.3,
                      loc="upper left", ncols=1)

    # Leyenda especial para single (posición + masa)
    if "single" in group_keys:
        ax_s = axes[group_keys.index("single")]
        from matplotlib.lines import Line2D
        positions = sorted({d["meta"]["pos"] for d in groups_data["single"] if d.get("meta")})
        pos_handles = [Line2D([0],[0], color=cm.get_cmap("cool")(i/max(len(positions)-1,1)),
                               lw=2, label=f"{p:.0f} AU")
                       for i, p in enumerate(positions)]
        mass_handles = [
            Line2D([0],[0], color="gray", ls=":", lw=1.5, label="Super-Earth"),
            Line2D([0],[0], color="gray", ls="-.", lw=1.5, label="Neptune"),
            Line2D([0],[0], color="gray", ls="--", lw=1.5, label="Saturn"),
            Line2D([0],[0], color="gray", ls="-",  lw=1.5, label="Jupiter"),
        ]
        l1 = ax_s.legend(handles=pos_handles, title="Pos. gap", fontsize=7,
                          loc="upper left", framealpha=0.35)
        ax_s.add_artist(l1)
        ax_s.legend(handles=mass_handles, title="Masa gap", fontsize=7,
                     loc="lower right", framealpha=0.35)

    fig.suptitle(f"Evolución de masa  —  embrión @ {r_planet_au:.1f} AU",
                 fontsize=14, fontweight="bold", y=1.01)
    fig.tight_layout()
    out = os.path.join(savedir, "mass_grouped.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"  [OK] {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# Figura 4: Flujo de pebbles en la snowline
# ═══════════════════════════════════════════════════════════════════════════════

def plot_pebble_flux_snow(all_data, savedir):
    """Calcula Mdot_pebble en r_snow para cada snapshot y run."""
    fig, ax = plt.subplots(figsize=(11, 5))

    group_color = {"single": "#4fc3f7", "multi": "#ff8a65", "sinusoidal": "#a5d6a7"}
    cmap_s = cm.get_cmap("cool")

    for d in all_data:
        pa3 = d.get("pa3")
        if pa3 is None:
            continue
        rsnow_arr = pa3.rsnow.get("H2O", np.array([]))
        if len(rsnow_arr) == 0:
            continue

        Nt = pa3.Nt
        t_yr  = pa3.times / YR
        flux  = np.zeros(Nt)

        for it in range(Nt):
            rsnow = rsnow_arr[it]
            if not np.isfinite(rsnow) or rsnow <= 0:
                continue
            # Sigma_dust total en la snowline
            sig_d = np.interp(np.log(rsnow), np.log(pa3.r),
                              pa3.dust["Sigma_total"][it])
            # v_r del polvo (bin grande = índice -1)
            vr    = np.interp(np.log(rsnow), np.log(pa3.r),
                              pa3.dust["vr"][it, :, -1])
            flux[it] = 2 * np.pi * rsnow * sig_d * abs(vr)  # g/s

        mask = flux > 0
        if mask.sum() < 2:
            continue

        color = group_color.get(d["group"], "#aaa")
        ax.plot(t_yr[mask], flux[mask] / (M_EARTH / YR),
                color=color, lw=1.2, alpha=0.55)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Tiempo [años]")
    ax.set_ylabel(r"$\dot{M}_{\rm peb}$ en snowline [$M_\oplus$ / yr]")
    ax.set_title("Flujo de pebbles en la snowline de H₂O  —  todas las runs", pad=10)

    from matplotlib.patches import Patch
    legend_els = [Patch(facecolor=v, label=k, alpha=0.8)
                  for k, v in group_color.items()]
    ax.legend(handles=legend_els, fontsize=9, framealpha=0.35)

    fig.tight_layout()
    out = os.path.join(savedir, "pebble_flux_snow.png")
    fig.savefig(out, bbox_inches="tight")
    print(f"  [OK] {out}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════════════
# Pipeline principal
# ═══════════════════════════════════════════════════════════════════════════════

def run_pipeline(dataroot, savedir, r_planet_au, t_min_yr):
    os.makedirs(savedir, exist_ok=True)
    runs = collect_runs(dataroot)
    print(f"\n  {len(runs)} runs encontradas en {dataroot}\n")

    single_runs, multi_runs, sinus_runs = group_runs(runs)

    all_data = []
    for group_list, gkey in [(single_runs, "single"),
                              (multi_runs,  "multi"),
                              (sinus_runs,  "sinusoidal")]:
        for d in group_list:
            print(f"  PA3 -> {d['name']} ...")
            try:
                pa3, hist = run_pa3(d["path"], r_planet_au, t_min_yr)
                d["hist"]  = hist
                d["pa3"]   = pa3
                d["group"] = gkey
                all_data.append(d)
                pa3.summary({r_planet_au: hist})
            except Exception as e:
                print(f"  [ERR] {d['name']}: {e}")

    print(f"\n  Generando figuras en {savedir}/ ...")

    # Fig 1: heatmap (solo single)
    single_done = [d for d in all_data if d["group"] == "single"]
    if single_done:
        plot_heatmap(single_done, savedir, r_planet_au)

    # Fig 2: bar chart ordenado
    plot_bar_fh2o(all_data, savedir, r_planet_au)

    # Fig 3: evolución agrupada
    groups_data = {
        "single":     [d for d in all_data if d["group"] == "single"],
        "multi":      [d for d in all_data if d["group"] == "multi"],
        "sinusoidal": [d for d in all_data if d["group"] == "sinusoidal"],
    }
    plot_mass_grouped(groups_data, savedir, r_planet_au)

    # Fig 4: flujo en la snowline
    plot_pebble_flux_snow(all_data, savedir)

    print(f"\n  Listo. Figuras guardadas en {savedir}/")


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Análisis PA3 para data_1myr/")
    p.add_argument("--dataroot",  default="data/1myr")
    p.add_argument("--savedir",   default="figures/1myr")
    p.add_argument("--r_planet",  type=float, default=2.0)
    p.add_argument("--t_min_yr",  type=float, default=100.0)
    args = p.parse_args()

    run_pipeline(args.dataroot, args.savedir, args.r_planet, args.t_min_yr)

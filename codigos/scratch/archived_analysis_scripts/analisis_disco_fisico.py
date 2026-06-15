# -*- coding: utf-8 -*-
"""
analisis_disco_fisico.py
------------------------
Análisis de parámetros físicos del disco extraídos directamente de los HDF5 de TripodPy.

Genera cuatro familias de gráficos:
  1. eta(r)          — gradiente de presión vs. distancia, para distintas M_gap y r_gap
  2. a_max(r)        — tamaño máximo de grano vs. distancia con isolíneas de Stokes fijo
  3. Flujo de pebbles (Fi.tot) alrededor del gap, comparativa de masas
  4. Mapa de calor: eta_max alrededor del gap vs. (M_gap, r_gap)

Fuente de datos: 10Myr_round0.1 (alpha = 0.001 fijo, múltiples r_gap y M_gap)

Autor: generado con Antigravity
"""

import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm, Normalize, BoundaryNorm
from matplotlib.cm import ScalarMappable
import h5py

# ===========================================================================
# CONSTANTES CGS
# ===========================================================================
AU_cm   = 1.495978707e13   # 1 AU en cm
M_jup_g = 1.898e30         # 1 M_Jup en gramos
year_s  = 3.15576e7        # 1 año en segundos

# ===========================================================================
# CONFIGURACIÓN
# ===========================================================================
BASE_DIR    = r"data/runs/10Myr_round0.1"
FIG_DIR     = r"data/figures/disco_fisico"
CACHE_FILE  = r"data/runs/disco_fisico_cache.pkl"
SNAPSHOT    = -1            # -1 = último snapshot (data0099.hdf5)

# Isolíneas de Stokes a trazar en el gráfico de tamaño de grano
STOKES_LINES = [0.001, 0.01, 0.05, 0.1, 0.3, 1.0]

# Radio de la ventana alrededor del gap para extraer eta_peak y flujo
WINDOW_FAC  = 2.0           # ventana = r_gap * [1/WINDOW_FAC , WINDOW_FAC]

os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    'font.size': 13,
    'axes.labelsize': 13,
    'axes.titlesize': 13,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'font.family': 'sans-serif',
})

# ===========================================================================
# UTILIDADES DE I/O
# ===========================================================================

def parse_run_name(run_name):
    """Extrae r_gap [AU], M_gap [M_Jup], alpha del nombre del directorio."""
    parts = run_name.split('_')
    r_val = float(parts[1][1:])
    m_val = float(parts[2][1:])
    a_val = float(parts[3][1:])
    return r_val, m_val, a_val


def get_completed_runs(base_dir=BASE_DIR):
    runs = sorted(glob.glob(os.path.join(base_dir, "run_*")))
    return [r for r in runs if os.path.exists(os.path.join(r, "data0099.hdf5"))]


def load_snapshot(run_path, snap_idx=SNAPSHOT):
    """
    Carga un snapshot HDF5 y devuelve un dict con las cantidades físicas relevantes.
    snap_idx: índice del snapshot (0-based). -1 = último disponible.
    """
    files = sorted(glob.glob(os.path.join(run_path, "data*.hdf5")))
    if not files:
        raise FileNotFoundError(f"No hay HDF5 en {run_path}")
    fpath = files[snap_idx]

    with h5py.File(fpath, "r") as f:
        r       = f["grid/r"][:]           # [cm]
        ri      = f["grid/ri"][:]          # interfaces [cm]
        OmegaK  = f["grid/OmegaK"][:]     # [1/s]
        t       = float(f["t"][()])        # [s]

        # GAS
        Sigma_g = f["gas/Sigma"][:]        # [g/cm²]
        eta     = f["gas/eta"][:]          # dimensionless
        cs      = f["gas/cs"][:]           # [cm/s]
        Hp      = f["gas/Hp"][:]           # [cm]
        P       = f["gas/P"][:]            # [g/cm/s²]
        rho_g   = f["gas/rho"][:]          # [g/cm³]
        alpha_g = f["gas/alpha"][:]        # dimensionless
        v_visc  = f["gas/v/visc"][:]       # [cm/s]
        v_rad_g = f["gas/v/rad"][:]        # [cm/s]
        Fi_gas  = f["gas/Fi"][:]           # [g/cm/s]

        # DUST
        a_all   = f["dust/a"][:]           # [cm] shape (Nr, 5)
        a_max   = f["dust/s/max"][:]       # [cm] shape (Nr,)
        St_all  = f["dust/St"][:]          # shape (Nr, 5)
        Sigma_d = f["dust/Sigma"][:]       # shape (Nr, 2)
        eps     = f["dust/eps"][:]         # dust-to-gas ratio
        rhos    = f["dust/rhos"][:, 0]     # bulk density [g/cm³]
        Fi_tot  = f["dust/Fi/tot"][:]      # shape (Nr+1, 2)
        v_frag  = f["dust/v/frag"][:]      # [cm/s]

    r_au   = r / AU_cm
    ri_au  = ri / AU_cm
    t_yr   = t / year_s

    # Flujo de pebbles total: suma sobre bins de masa, en los centros de celda
    # Fi_tot tiene shape (Nr+1, 2) en interfaces; promediamos a centros
    Fi_pebble_if = Fi_tot[:, 1]   # bin grande (pebbles)
    # Interpolamos a centros: promedio entre interfaz izq y der
    Fi_pebble = 0.5 * (Fi_pebble_if[:-1] + Fi_pebble_if[1:])

    # Stokes del grano mayor (índice 4 = a_max)
    St_max = St_all[:, 4]

    return {
        "r_au":       r_au,
        "ri_au":      ri_au,
        "t_yr":       t_yr,
        "eta":        eta,
        "cs":         cs,
        "Hp":         Hp,
        "P":          P,
        "OmegaK":     OmegaK,
        "Sigma_g":    Sigma_g,
        "rho_g":      rho_g,
        "alpha":      alpha_g,
        "v_visc":     v_visc,
        "v_rad_g":    v_rad_g,
        "Fi_gas":     Fi_gas,
        "a_max":      a_max,
        "St_max":     St_max,
        "Sigma_d":    Sigma_d,
        "eps":        eps,
        "rhos":       rhos,
        "Fi_pebble":  Fi_pebble,
        "v_frag":     v_frag,
    }

# ===========================================================================
# EXTRACCIÓN Y CACHE
# ===========================================================================

def extract_all(base_dir=BASE_DIR, cache_file=CACHE_FILE, snap_idx=SNAPSHOT, force=False):
    if os.path.exists(cache_file) and not force:
        print(f"Cargando cache: {cache_file}")
        with open(cache_file, "rb") as f:
            return pickle.load(f)

    runs = get_completed_runs(base_dir)
    print(f"Procesando {len(runs)} simulaciones completas en {base_dir}...")
    data = []

    for i, rpath in enumerate(runs):
        run_name = os.path.basename(rpath)
        r_gap, M_gap, alpha = parse_run_name(run_name)
        try:
            snap = load_snapshot(rpath, snap_idx)
        except Exception as e:
            print(f"  [{i+1}] ERROR en {run_name}: {e}")
            continue

        # Índice del gap en la grilla (celda más cercana a r_gap)
        i_gap = np.argmin(np.abs(snap["r_au"] - r_gap))

        # Ventana de ±WINDOW_FAC alrededor del gap
        r_in  = r_gap / WINDOW_FAC
        r_out = r_gap * WINDOW_FAC
        mask_win = (snap["r_au"] >= r_in) & (snap["r_au"] <= r_out)

        # eta_peak en la zona del gap (máximo del gradiente, presión bump exterior)
        eta_win    = snap["eta"][mask_win]
        eta_peak   = np.max(eta_win) if len(eta_win) > 0 else np.nan

        # Flujo de pebbles en el borde exterior del gap
        # Tomamos la celda justo en r_gap (flujo entrante al gap)
        Fi_gap_ext = snap["Fi_pebble"][i_gap] if i_gap < len(snap["Fi_pebble"]) else np.nan

        # a_max y St: tomamos el valor en el pressure bump exterior al gap
        # (celda con eta máximo, que es donde los granos se acumulan)
        if len(eta_win) > 0:
            i_bump = np.argmax(snap["eta"][mask_win])
            r_win  = snap["r_au"][mask_win]
            i_bump_global = np.argmin(np.abs(snap["r_au"] - r_win[i_bump]))
        else:
            i_bump_global = i_gap
        a_gap_cm   = snap["a_max"][i_bump_global]
        St_gap     = snap["St_max"][i_bump_global]

        data.append({
            "run_name":   run_name,
            "r_gap":      r_gap,
            "M_gap":      M_gap,
            "alpha":      alpha,
            "t_yr":       snap["t_yr"],
            # Perfiles completos
            "r_au":       snap["r_au"],
            "eta":        snap["eta"],
            "a_max":      snap["a_max"],
            "St_max":     snap["St_max"],
            "Sigma_g":    snap["Sigma_g"],
            "Sigma_d":    snap["Sigma_d"],
            "eps":        snap["eps"],
            "Fi_pebble":  snap["Fi_pebble"],
            "rhos":       snap["rhos"],
            "cs":         snap["cs"],
            "OmegaK":     snap["OmegaK"],
            # Escalares útiles
            "eta_peak":   eta_peak,
            "Fi_gap_ext": Fi_gap_ext,
            "a_gap_cm":   a_gap_cm,
            "St_gap":     St_gap,
        })
        print(f"  [{i+1}/{len(runs)}] {run_name} OK  "
              f"| eta_peak={eta_peak:.4f}  a_max_gap={a_gap_cm*100:.3f} mm  St={St_gap:.4f}")

    with open(cache_file, "wb") as f:
        pickle.dump(data, f)
    print(f"Cache guardado en {cache_file}")
    return data


# ===========================================================================
# PALETA DIVERGENTE CENTRADA EN M_peak
# ===========================================================================
import colorsys

def make_mgap_palette(all_M_gaps, M_peak=0.1):
    """
    Devuelve un dict {M_gap: color} donde M_gap == M_peak es el más
    saturado/oscuro, y los colores se atenúan simétricamente (en log)
    alejándose del peak. Usa azules para M ≤ M_peak, rojos para M > M_peak.
    """
    M_peak_log = np.log10(M_peak)
    log_dists   = {m: abs(np.log10(m) - M_peak_log) for m in all_M_gaps}
    max_dist    = max(log_dists.values()) if log_dists else 1.0

    palette = {}
    for m in all_M_gaps:
        nd = log_dists[m] / max_dist        # 0 = en el peak, 1 = más lejos
        if m <= M_peak:
            hue = 0.62                      # azul
        else:
            hue = 0.04                      # rojo-anaranjado
        lightness  = 0.18 + 0.58 * nd      # oscuro en peak, claro lejos
        saturation = 0.90 - 0.30 * nd      # muy saturado en peak
        r, g, b = colorsys.hls_to_rgb(hue, lightness, saturation)
        palette[m] = (r, g, b)
    return palette


def add_mgap_legend(ax, palette, all_M_gaps, M_peak=0.1):
    """
    Agrega al eje ax una leyenda compacta de colores de M_gap
    en el lado derecho sin obstruir los datos.
    """
    from matplotlib.lines import Line2D
    handles = []
    labels  = []
    for m in sorted(all_M_gaps):
        lw   = 2.5 if m == M_peak else 1.8
        ls   = "-"
        handles.append(Line2D([0], [0], color=palette[m], lw=lw, ls=ls))
        label = rf"${m}\,M_{{\rm Jup}}$"
        if m == M_peak:
            label += r" $\star$"
        labels.append(label)
    leg = ax.legend(handles, labels, loc="upper right", fontsize=8,
                    framealpha=0.7, edgecolor="gray",
                    title=r"$M_{\rm gap}$", title_fontsize=8,
                    handlelength=1.5, labelspacing=0.3)
    return leg


# ===========================================================================
# PLOT 1: eta(r) agrupado por r_gap, coloreado por M_gap
# ===========================================================================

def plot_eta_vs_r(data, alpha_filter=0.001, r_gaps_show=None, M_peak=0.1):
    """
    Perfiles de eta(r) coloreados con paleta divergente centrada en M_peak.
    Un panel por r_gap. Labels de ejes solo donde corresponde.
    Leyenda de M_gap dentro del plot en vez de colorbar externo.
    """
    subset = [d for d in data if abs(d["alpha"] - alpha_filter) < 1e-6]
    if not subset:
        print(f"No hay datos para alpha={alpha_filter}")
        return

    all_r_gaps = sorted(set(d["r_gap"] for d in subset))
    if r_gaps_show is not None:
        all_r_gaps = [r for r in all_r_gaps if r in r_gaps_show]

    all_M_gaps = sorted(set(d["M_gap"] for d in subset))
    palette    = make_mgap_palette(all_M_gaps, M_peak)

    cols = 3
    rows = (len(all_r_gaps) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows),
                             constrained_layout=True)
    fig.suptitle(
        rf"Gradiente de presión $\eta(r)$ — $\alpha = {alpha_filter}$"
        rf" — $t = {subset[0]['t_yr']/1e6:.1f}$ Myr",
        fontsize=16
    )
    axes = np.array(axes).flatten()

    for ax_i, rg in enumerate(all_r_gaps):
        ax   = axes[ax_i]
        row  = ax_i // cols
        col  = ax_i % cols
        runs = sorted([d for d in subset if d["r_gap"] == rg],
                      key=lambda x: x["M_gap"])

        for d in runs:
            color = palette[d["M_gap"]]
            lw    = 2.2 if d["M_gap"] == M_peak else 1.6
            ax.plot(d["r_au"], d["eta"], color=color, lw=lw, alpha=0.90)

        ax.axvline(rg, color="gray", ls="--", lw=0.9, alpha=0.55)
        ax.set_xlim([rg / WINDOW_FAC * 0.6, rg * WINDOW_FAC * 1.4])

        all_eta_win = np.concatenate([
            d["eta"][(d["r_au"] >= rg / WINDOW_FAC * 0.6)
                     & (d["r_au"] <= rg * WINDOW_FAC * 1.4)]
            for d in runs
        ])
        if len(all_eta_win):
            ax.set_ylim([0, np.nanpercentile(all_eta_win, 98) * 1.25])

        ax.set_title(f"Gap @ {rg} AU", fontsize=12)
        ax.grid(True, alpha=0.22, lw=0.5)

        # Labels de eje: Y solo columna izquierda, X solo última fila
        if col == 0:
            ax.set_ylabel(r"$\eta$")
        else:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)

        # Calcular última fila activa de esta columna
        last_active_row = max(
            (i // cols) for i in range(len(all_r_gaps)) if i % cols == col
        )
        if row == last_active_row:
            ax.set_xlabel(r"$r$ [AU]")
        else:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)

        # Leyenda de colores solo en el primer panel (top-left)
        if ax_i == 0:
            add_mgap_legend(ax, palette, all_M_gaps, M_peak)

    for j in range(len(all_r_gaps), len(axes)):
        axes[j].set_visible(False)

    out = os.path.join(FIG_DIR, f"eta_vs_r_alpha{alpha_filter}.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Guardado: {out}")



# ===========================================================================
# PLOT 2: a_max(r) + isolíneas de St fijo
# ===========================================================================

def stokes_to_amax(St_target, Sigma_g, rhos):
    """
    Tamaño de grano en régimen Epstein para un St dado:
      a = (2/pi) * St * Sigma_g / rhos
    """
    return (2.0 / np.pi) * St_target * Sigma_g / rhos


def plot_amax_vs_r_with_stokes(data, alpha_filter=0.001, r_gaps_show=None,
                                stokes_lines=STOKES_LINES, M_peak=0.1):
    """
    a_max(r) con paleta divergente centrada en M_peak.
    Isolíneas de St en negro semitransparente, etiquetas al lado derecho.
    Labels de ejes solo donde corresponde. Leyenda de M_gap inline.
    """
    subset = [d for d in data if abs(d["alpha"] - alpha_filter) < 1e-6]
    if not subset:
        return

    all_r_gaps = sorted(set(d["r_gap"] for d in subset))
    if r_gaps_show is not None:
        all_r_gaps = [r for r in all_r_gaps if r in r_gaps_show]

    all_M_gaps = sorted(set(d["M_gap"] for d in subset))
    palette    = make_mgap_palette(all_M_gaps, M_peak)

    cols = 3
    rows = (len(all_r_gaps) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows),
                             constrained_layout=True)
    fig.suptitle(
        rf"Tamaño máximo de grano $a_{{\rm max}}(r)$ — $\alpha = {alpha_filter}$",
        fontsize=16
    )
    axes = np.array(axes).flatten()

    for ax_i, rg in enumerate(all_r_gaps):
        ax   = axes[ax_i]
        row  = ax_i // cols
        col  = ax_i % cols
        runs = sorted([d for d in subset if d["r_gap"] == rg],
                      key=lambda x: x["M_gap"])

        # --- Perfiles a_max ---
        for d in runs:
            color  = palette[d["M_gap"]]
            lw     = 2.2 if d["M_gap"] == M_peak else 1.6
            a_plot = np.where(d["a_max"] > 1e-6, d["a_max"], np.nan)
            ax.semilogy(d["r_au"], a_plot, color=color, lw=lw, alpha=0.90)

        # --- Isolíneas de Stokes (negro, alpha=0.3) ---
        # Usamos el Sigma_g del run de menor M_gap (sin gap = unperturbed)
        ref    = runs[0]
        r_plot = ref["r_au"]
        x_max  = rg * WINDOW_FAC * 1.4   # borde derecho visible

        for st in stokes_lines:
            a_iso = stokes_to_amax(st, ref["Sigma_g"], ref["rhos"])
            ax.semilogy(r_plot, a_iso, color="black", lw=1.0,
                        ls="--", alpha=0.3, zorder=2)
            # Etiqueta al borde derecho del panel
            idx_r = np.argmin(np.abs(r_plot - x_max * 0.97))
            y_val = a_iso[idx_r]
            ylim  = ax.get_ylim()
            # Solo etiquetar si el valor está dentro del rango Y visible
            # (se aplica después de set_ylim; guardamos para aplicar post-hoc)
            ax.text(x_max * 0.98, y_val, f" St={st}",
                    fontsize=7.5, color="black", alpha=0.55,
                    ha="right", va="center",
                    clip_on=True)

        ax.axvline(rg, color="gray", ls="--", lw=0.9, alpha=0.55)
        ax.set_xlim([rg / WINDOW_FAC * 0.6, x_max])

        # Rango Y físico
        all_a_win = np.concatenate([
            d["a_max"][(d["r_au"] >= rg / WINDOW_FAC * 0.6)
                       & (d["r_au"] <= x_max)]
            for d in runs
        ])
        all_a_win = all_a_win[all_a_win > 1e-6]
        if len(all_a_win) > 0:
            ax.set_ylim([np.nanmin(all_a_win) * 0.4, np.nanmax(all_a_win) * 4.0])

        ax.set_title(f"Gap @ {rg} AU", fontsize=12)
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(
            lambda x, _: f"{x*10:.2g} mm" if x < 0.1 else f"{x:.2g} cm"
        ))
        ax.grid(True, alpha=0.22, lw=0.5, which="both")

        # Labels de eje: Y solo columna izquierda, X solo última fila activa
        if col == 0:
            ax.set_ylabel(r"$a_{\rm max}$")
        else:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)

        last_active_row = max(
            (i // cols) for i in range(len(all_r_gaps)) if i % cols == col
        )
        if row == last_active_row:
            ax.set_xlabel(r"$r$ [AU]")
        else:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)

        # Leyenda de M_gap solo en primer panel
        if ax_i == 0:
            add_mgap_legend(ax, palette, all_M_gaps, M_peak)

    for j in range(len(all_r_gaps), len(axes)):
        axes[j].set_visible(False)

    out = os.path.join(FIG_DIR, f"amax_stokes_alpha{alpha_filter}.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Guardado: {out}")


# ===========================================================================
# PLOT 3: Flujo de pebbles alrededor del gap
# ===========================================================================

def plot_pebble_flux_around_gap(data, alpha_filter=0.001, r_gaps_show=None,
                                M_peak=0.1):
    """
    |Fi_pebble|(r) normalizado por el run de M_gap mínimo (casi sin gap),
    con paleta divergente centrada en M_peak. Labels de ejes solo donde
    corresponde. Leyenda inline en vez de colorbar.
    """
    subset = [d for d in data if abs(d["alpha"] - alpha_filter) < 1e-6]
    if not subset:
        return

    all_r_gaps = sorted(set(d["r_gap"] for d in subset))
    if r_gaps_show is not None:
        all_r_gaps = [r for r in all_r_gaps if r in r_gaps_show]

    all_M_gaps = sorted(set(d["M_gap"] for d in subset))
    palette    = make_mgap_palette(all_M_gaps, M_peak)

    cols = 3
    rows = (len(all_r_gaps) + cols - 1) // cols
    fig, axes = plt.subplots(rows, cols, figsize=(5.5 * cols, 4.5 * rows),
                             constrained_layout=True)
    fig.suptitle(
        rf"Flujo de pebbles $|\dot{{F}}_{{pebbles}}|(r)$ — $\alpha = {alpha_filter}$",
        fontsize=16
    )
    axes = np.array(axes).flatten()

    for ax_i, rg in enumerate(all_r_gaps):
        ax   = axes[ax_i]
        row  = ax_i // cols
        col  = ax_i % cols
        runs = sorted([d for d in subset if d["r_gap"] == rg],
                      key=lambda x: x["M_gap"])

        # Referencia: run de menor M_gap (casi sin gap)
        Fi_ref      = np.abs(runs[0]["Fi_pebble"])
        Fi_ref_safe = np.where(Fi_ref > 0, Fi_ref, np.nan)

        for d in runs:
            color    = palette[d["M_gap"]]
            lw       = 2.2 if d["M_gap"] == M_peak else 1.6
            r_au     = d["r_au"]
            Fi       = np.abs(d["Fi_pebble"])
            Fi_norm  = Fi / Fi_ref_safe
            mask_win = (r_au >= rg / WINDOW_FAC * 0.5) & (r_au <= rg * WINDOW_FAC * 1.5)
            ax.semilogy(r_au[mask_win], Fi_norm[mask_win],
                        color=color, lw=lw, alpha=0.90)

        ax.axvline(rg, color="gray", ls="--", lw=0.9, alpha=0.55)
        ax.axhline(1.0, color="k",    ls=":",  lw=1.0, alpha=0.4)
        ax.set_xlim([rg / WINDOW_FAC * 0.5, rg * WINDOW_FAC * 1.5])
        ax.set_title(f"Gap @ {rg} AU", fontsize=12)
        ax.grid(True, alpha=0.22, lw=0.5, which="both")

        # Labels de eje: Y solo columna izquierda, X solo última fila activa
        if col == 0:
            ax.set_ylabel(r"$|\dot{F}_{\rm pebble}| / |\dot{F}_{\rm ref}|$")
        else:
            ax.set_ylabel("")
            ax.tick_params(labelleft=False)

        last_active_row = max(
            (i // cols) for i in range(len(all_r_gaps)) if i % cols == col
        )
        if row == last_active_row:
            ax.set_xlabel(r"$r$ [AU]")
        else:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)

        # Leyenda de M_gap solo en primer panel
        if ax_i == 0:
            add_mgap_legend(ax, palette, all_M_gaps, M_peak)

    for j in range(len(all_r_gaps), len(axes)):
        axes[j].set_visible(False)

    out = os.path.join(FIG_DIR, f"pebble_flux_alpha{alpha_filter}.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Guardado: {out}")



# ===========================================================================
# PLOT 4: Mapa de calor eta_peak vs (r_gap, M_gap) — el gráfico estrella
# ===========================================================================

def plot_eta_peak_heatmap(data, alpha_filter=0.001):
    """
    Heatmap 2D: eje X = r_gap [AU], eje Y = M_gap [M_Jup],
    color = eta_peak (el máximo de eta en la ventana alrededor del gap).
    """
    subset = [d for d in data if abs(d["alpha"] - alpha_filter) < 1e-6]
    if not subset:
        return

    r_gaps = sorted(set(d["r_gap"] for d in subset))
    M_gaps = sorted(set(d["M_gap"] for d in subset))

    Z = np.full((len(M_gaps), len(r_gaps)), np.nan)
    for d in subset:
        j = r_gaps.index(d["r_gap"])
        i = M_gaps.index(d["M_gap"])
        Z[i, j] = d["eta_peak"]

    fig, ax = plt.subplots(figsize=(10, 6))
    vmin, vmax = np.nanmin(Z), np.nanmax(Z)
    im = ax.imshow(Z, origin="lower", aspect="auto",
                   cmap="inferno",
                   vmin=vmin, vmax=vmax,
                   extent=[-0.5, len(r_gaps) - 0.5, -0.5, len(M_gaps) - 0.5])

    # Texto en cada celda
    for i in range(len(M_gaps)):
        for j in range(len(r_gaps)):
            val = Z[i, j]
            if not np.isnan(val):
                txt_color = "white" if val < (vmin + 0.6 * (vmax - vmin)) else "black"
                ax.text(j, i, f"{val:.4f}", ha="center", va="center",
                        color=txt_color, fontsize=9)

    ax.set_xticks(range(len(r_gaps)))
    ax.set_xticklabels([str(r) for r in r_gaps])
    ax.set_yticks(range(len(M_gaps)))
    ax.set_yticklabels([str(m) for m in M_gaps])
    ax.set_xlabel(r"$r_{\rm gap}$ [AU]")
    ax.set_ylabel(r"$M_{\rm gap}$ [$M_{\rm Jup}$]")
    ax.set_title(
        rf"$\eta_{{\rm peak}}$ alrededor del gap — $\alpha = {alpha_filter}$  "
        rf"($t \approx {subset[0]['t_yr']/1e6:.1f}$ Myr)"
    )

    cbar = fig.colorbar(im, ax=ax, pad=0.02)
    cbar.set_label(r"$\eta_{\rm peak}$")

    # Línea horizontal en M_gap ≈ 0.1 M_Jup (transición conocida)
    if 0.1 in M_gaps:
        i_trans = M_gaps.index(0.1)
        ax.axhline(i_trans - 0.5, color="cyan", ls="--", lw=2.0, alpha=0.7,
                   label=r"Transición $M_{\rm gap}=0.1\,M_{\rm Jup}$")
        ax.legend(loc="upper right", fontsize=10)

    plt.tight_layout()
    out = os.path.join(FIG_DIR, f"heatmap_eta_peak_alpha{alpha_filter}.png")
    plt.savefig(out, dpi=180, bbox_inches="tight")
    plt.close()
    print(f"  Guardado: {out}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    print("=" * 60)
    print("ANÁLISIS FÍSICO DEL DISCO (TripodPy HDF5)")
    print("=" * 60)

    data = extract_all(BASE_DIR, CACHE_FILE, SNAPSHOT)

    if not data:
        print("No se encontraron datos. Revisá BASE_DIR y los archivos HDF5.")
        return

    print(f"\nTotal de simulaciones procesadas: {len(data)}")
    alphas = sorted(set(d["alpha"] for d in data))
    print(f"Alphas disponibles: {alphas}")

    for alpha in alphas:
        n = sum(1 for d in data if abs(d["alpha"] - alpha) < 1e-6)
        print(f"\n--- alpha = {alpha}  ({n} corridas) ---")
        print("  Generando gráficos...")

        plot_eta_vs_r(data, alpha)
        plot_amax_vs_r_with_stokes(data, alpha)
        plot_pebble_flux_around_gap(data, alpha)
        plot_eta_peak_heatmap(data, alpha)

    print(f"\n¡Todos los gráficos guardados en: {FIG_DIR}!")


if __name__ == "__main__":
    main()

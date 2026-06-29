# -*- coding: utf-8 -*-
"""
plot_figura_vfrag_alpha.py
==========================
Genera tres figuras para §4.2 "Crecimiento de un Embrión Individual":

  Figura A1 — Gap único (1×3):
      Columnas: v_frag = 1, 3, 10 m/s
      Colorbar (inferno): M_gap  [M_Jup]

  Figura A2 — Sinusoidal (1×3):
      Columnas: v_frag = 1, 3, 10 m/s
      Colorbar (viridis): Amplitud A

  Figura B — Regímenes de M_gap (2×2):
      Paneles: M_gap = 0.01, 0.1, 0.3, 0.5
      Colorbar (RdYlBu_r): alpha

  Ejes comunes:
      X: log-scale, labels log₁₀(t [Myr]) → -3, -2, -1, 0, 1
      Y: log-scale, labels log₁₀(Mₑₘᵦ [M⊕])
      Título: alpha = 0.001, M_{emb,0} = 10^{-4} M_earth

Salida:
    data/figures/paper/figura_A1_gap_unico.png
    data/figures/paper/figura_A2_sinusoidal.png
    data/figures/paper/figura_B_regimenes_mgap.png
"""

import os, sys, io, pickle, contextlib, re, math
import numpy as np

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
else:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm, Normalize
from matplotlib.cm import ScalarMappable

# ── PA3Py ─────────────────────────────────────────────────────────────────────
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
try:
    from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
    import dustpy.constants as c
    M_EARTH = c.M_earth
    YEAR    = c.year
except ImportError as e:
    print(f"[AVISO] {e}")
    M_EARTH = 5.972e27; YEAR = 3.156e7
    PebbleAccretionModule3 = None

# ── Parámetros ────────────────────────────────────────────────────────────────
M0_EARTH        = 1e-4
ALPHA_CANONICAL = 0.001
VFRAG_CANONICAL = 10
RGAP_FIXED      = 10.0
NGAP_SIN        = 10
ALPHAS_B        = [0.0001, 0.0005, 0.001, 0.003]
MGAPS_B         = [0.01, 0.1, 0.3, 0.5]
VFRAGS_A        = [1, 3, 10]
BASE_RUNS       = 'data/runs'
OUT_DIR         = 'data/figures/paper'
os.makedirs(OUT_DIR, exist_ok=True)

STYLES = ['-', '--', '-.', ':', (0,(5,5)), (0,(3,5,1,5))]

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 13,
    'axes.labelsize': 15,
    'xtick.labelsize': 13,
    'ytick.labelsize': 13,
    'legend.fontsize': 11,
})

# Formateadores de ejes (log-scale pero labels muestran log10)
def _log_fmt_x(x, pos):
    if x <= 0: return ""
    v = np.log10(x)
    if abs(v - round(v)) < 0.05:
        return f"${int(round(v))}$"
    return ""

def _log_fmt_y(x, pos):
    if x <= 0: return ""
    v = np.log10(x)
    if abs(v - round(v)) < 0.05:
        return f"${int(round(v))}$"
    return ""

XFMT = ticker.FuncFormatter(_log_fmt_x)
YFMT = ticker.FuncFormatter(_log_fmt_y)

XLABEL = r"$\log_{10}(t\,[\mathrm{Myr}])$"
YLABEL = r"$\log_{10}(M_{\mathrm{emb}}\,[M_\oplus])$"
TITLE_ALPHA = (rf"$\alpha = {ALPHA_CANONICAL}$,"
               r"  $M_{\mathrm{emb},0} = 10^{-4}\,M_\oplus$")


# ══════════════════════════════════════════════════════════════════════════════
# EXTRACCIÓN
# ══════════════════════════════════════════════════════════════════════════════

def _extract_run(rpath):
    if PebbleAccretionModule3 is None:
        raise RuntimeError("PA3Py no disponible.")
    with open(os.devnull, 'w') as fnull, contextlib.redirect_stdout(fnull):
        pa3 = PebbleAccretionModule3.from_datadir(rpath, M_star=1.0)
        res = pa3.run_growth([1.0], M0_g=M0_EARTH * M_EARTH)
    hist = res[1.0]
    if len(hist) == 0: return None
    times_yr = hist[:,0] / YEAR
    mass_e   = hist[:,1] / M_EARTH
    rsnow_au = hist[:,4]
    m_iso_e  = hist[:,5] / M_EARTH
    ci = np.where(rsnow_au < 1.0)[0]
    t_cross = times_yr[ci[0]] if len(ci) > 0 else None
    return dict(times_yr=times_yr, mass_e=mass_e, m_iso_e=m_iso_e,
                t_cross_1au=t_cross)


def _load(cache_path, runs_dir, filter_fn, parse_fn):
    if os.path.exists(cache_path):
        print(f"  [CACHE] {cache_path}")
        with open(cache_path, 'rb') as f:
            return pickle.load(f)
    print(f"  [EXTRAYENDO] {runs_dir} ...")
    entries = []
    if not os.path.isdir(runs_dir): return entries
    for name in os.listdir(runs_dir):
        rpath = os.path.join(runs_dir, name)
        if not os.path.isdir(rpath): continue
        if not filter_fn(name): continue
        if not os.path.exists(os.path.join(rpath, "data0000.hdf5")): continue
        params = parse_fn(name)
        if params is None: continue
        try: hd = _extract_run(rpath)
        except Exception as e: print(f"    [ERR] {name}: {e}"); continue
        if hd is None: continue
        entries.append({**params, **hd, 'run_name': name})
    if entries:
        os.makedirs(os.path.dirname(cache_path), exist_ok=True)
        with open(cache_path, 'wb') as f: pickle.dump(entries, f)
        print(f"  [GUARDADO] {len(entries)} entradas")
    return entries


def load_general_r10(vfrag, alpha):
    runs_dir   = os.path.join(BASE_RUNS, f"vf_{vfrag}ms", "general")
    cache_path = os.path.join(runs_dir,
                              f"cache_general_r10_a_{alpha}_M0_{M0_EARTH}.pkl")
    def flt(n): return (f"_r{RGAP_FIXED}" in n or
                        f"_r{int(RGAP_FIXED)}." in n) and \
                       (n.endswith(f"_a{alpha}") or f"_a{alpha}_" in n)
    def parse(n):
        pts = n.split('_')
        try:
            r = float(pts[1][1:]); m = float(pts[2][1:]); a = float(pts[3][1:])
            if abs(r - RGAP_FIXED) > 0.5: return None
            return dict(r_gap=r, M_gap=m, alpha=a)
        except: return None
    return _load(cache_path, runs_dir, flt, parse)


def load_sinusoidal_ngap(vfrag, alpha, ngap):
    runs_dir   = os.path.join(BASE_RUNS, f"vf_{vfrag}ms", "sinusoidal")
    cache_path = os.path.join(runs_dir,
                              f"cache_sinusoidal_ngap{ngap}_a_{alpha}_M0_{M0_EARTH}.pkl")
    def flt(n): return (f"_ngap{ngap}_" in n) and \
                       (n.endswith(f"_a{alpha}") or f"_a{alpha}_" in n)
    def parse(n):
        try:
            ng = int(re.search(r'ngap(\d+)', n).group(1))
            av = float(re.search(r'_A([\d.]+)', n).group(1))
            al = float(re.search(r'_a([\d.]+)', n).group(1))
            return dict(n_gaps=ng, amp_val=av, alpha=al)
        except: return None
    return _load(cache_path, runs_dir, flt, parse)


# ══════════════════════════════════════════════════════════════════════════════
# UTILIDADES DE DIBUJO
# ══════════════════════════════════════════════════════════════════════════════

def _style_ax(ax, min_y=5e-5, max_y=12.0):
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlim([1e-3, 10.0]); ax.set_ylim([min_y, max_y])
    ax.grid(True, which='major', color='lightgray', ls='--', alpha=0.4)
    ax.grid(True, which='minor', color='whitesmoke', ls=':', alpha=0.25)
    ax.tick_params(direction='in', top=True, right=True,
                   which='major', length=6, width=1.2)
    ax.tick_params(direction='in', top=True, right=True,
                   which='minor', length=3, width=0.7)
    # X: solo ticks -2, -1, 0, 1  (evita superposicion con bordes de panel)
    ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
    ax.xaxis.set_major_formatter(XFMT)
    ax.yaxis.set_major_formatter(YFMT)
    # Umbral 0.1 M_earth
    ax.axhline(0.1, color='black', ls=':', alpha=0.85, lw=1.8, zorder=5)


def _draw_curves(ax, entries, color_fn, iso_color='saddlebrown', lw=3.0):
    t_cross_list = []
    for item in entries:
        color  = color_fn(item)
        t_myr  = item['times_yr'] / 1e6
        mass_e = item['mass_e']
        mask   = np.isfinite(t_myr) & np.isfinite(mass_e) & (mass_e > 0)
        if not mask.any(): continue
        tv = t_myr[mask]; mv = mass_e[mask]
        if tv[-1] < 10.0:
            tv = np.append(tv, 10.0); mv = np.append(mv, mv[-1])
        if 'm_iso_e' in item:
            ax.plot(t_myr[mask], item['m_iso_e'][mask],
                    color=iso_color, ls='--', alpha=0.5, lw=1.6, zorder=1)
        ax.plot(tv, mv, color=color, alpha=0.88, lw=lw, zorder=3)
        if item.get('t_cross_1au'):
            t_cross_list.append(item['t_cross_1au'] / 1e6)
    return t_cross_list


def _add_snowline(ax, tc_list):
    if tc_list:
        ax.axvline(np.median(tc_list), color='royalblue',
                   ls='--', alpha=0.75, lw=2.0, zorder=6)


def _add_colorbar(fig, axes_list, cmap, norm, label):
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    # Usa el ultimo ax para anclar la colorbar
    cbar = fig.colorbar(sm, ax=axes_list,
                        fraction=0.03, pad=0.02, aspect=30)
    cbar.set_label(label, fontsize=13)
    cbar.ax.tick_params(labelsize=11)
    return cbar


# ══════════════════════════════════════════════════════════════════════════════
# FIGURA A1 — Gap único 1×3
# ══════════════════════════════════════════════════════════════════════════════

def build_figura_A1():
    print("\n-- FIGURA A1: Gap unico 1x3 ------------------------------------")

    # Recopilar datos
    data = {vf: load_general_r10(vf, ALPHA_CANONICAL) for vf in VFRAGS_A}

    all_mgaps = sorted(set(e['M_gap'] for vf in VFRAGS_A for e in data[vf]))
    print(f"  M_gap encontrados: {all_mgaps}")

    cmap  = matplotlib.colormaps.get_cmap('inferno')
    norm  = LogNorm(vmin=max(min(all_mgaps), 1e-3), vmax=max(all_mgaps))
    def cfn(item): return cmap(norm(item['M_gap']))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5),
                             sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.02, right=0.88)

    for col, vf in enumerate(VFRAGS_A):
        ax = axes[col]
        _style_ax(ax)
        entries = data[vf]
        if entries:
            tc = _draw_curves(ax, sorted(entries, key=lambda x: x['M_gap']), cfn)
            _add_snowline(ax, tc)
        else:
            ax.text(0.5, 0.5, 'Sin datos', transform=ax.transAxes,
                    ha='center', va='center', color='gray')

        ax.set_title(rf"$v_{{\rm frag}} = {vf}$ m s$^{{-1}}$", fontsize=14, pad=6)
        ax.set_xlabel(XLABEL, fontsize=14)
        if col == 0:
            ax.set_ylabel(YLABEL, fontsize=14)

    # Leyenda M_iso y Umbral
    h = [Line2D([0],[0], color='saddlebrown', lw=1.6, ls='--', alpha=0.6),
         Line2D([0],[0], color='black', lw=1.8, ls=':', alpha=0.85),
         Line2D([0],[0], color='royalblue', lw=2.0, ls='--', alpha=0.75)]
    l = [r"$M_{\rm iso}$", r"$0.1\,M_\oplus$", "Snowline"]
    axes[0].legend(h, l, fontsize=10.5, loc='upper left', framealpha=0.85)

    # Colorbar M_gap
    _add_colorbar(fig, axes.tolist(), cmap, norm,
                  r"$M_{\rm gap}\,[M_{\rm Jup}]$")

    fig.suptitle(
        rf"$r_{{\mathrm{{gap}}}} = {int(RGAP_FIXED)}$ AU,  " + TITLE_ALPHA,
        fontsize=13, y=1.01
    )

    out_png = os.path.join(OUT_DIR, "figura_A1_gap_unico.png")
    
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    
    plt.close()
    print(f"  OK -> {out_png} ")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURA A2 — Sinusoidal 1×3
# ══════════════════════════════════════════════════════════════════════════════

def build_figura_A2():
    print("\n-- FIGURA A2: Sinusoidal 1x3 ------------------------------------")

    data = {vf: load_sinusoidal_ngap(vf, ALPHA_CANONICAL, NGAP_SIN)
            for vf in VFRAGS_A}

    all_amps = sorted(set(e['amp_val'] for vf in VFRAGS_A for e in data[vf]))
    print(f"  Amplitudes encontradas: {all_amps}")

    cmap = matplotlib.colormaps.get_cmap('viridis')
    norm = LogNorm(vmin=max(min(all_amps), 0.1), vmax=max(all_amps))
    def cfn(item): return cmap(norm(item['amp_val']))

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5),
                             sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.02, right=0.88)

    for col, vf in enumerate(VFRAGS_A):
        ax = axes[col]
        _style_ax(ax)
        entries = data[vf]
        if entries:
            tc = _draw_curves(ax, sorted(entries, key=lambda x: x['amp_val']), cfn)
            _add_snowline(ax, tc)
        else:
            ax.text(0.5, 0.5, 'Sin datos', transform=ax.transAxes,
                    ha='center', va='center', color='gray')

        ax.set_title(rf"$v_{{\rm frag}} = {vf}$ m s$^{{-1}}$", fontsize=14, pad=6)
        ax.set_xlabel(XLABEL, fontsize=14)
        if col == 0:
            ax.set_ylabel(YLABEL, fontsize=14)

    h = [Line2D([0],[0], color='saddlebrown', lw=1.6, ls='--', alpha=0.6),
         Line2D([0],[0], color='black', lw=1.8, ls=':', alpha=0.85),
         Line2D([0],[0], color='royalblue', lw=2.0, ls='--', alpha=0.75)]
    l = [r"$M_{\rm iso}$", r"Umbral $(0.1\,M_\oplus)$", "Snowline"]
    axes[0].legend(h, l, fontsize=10.5, loc='upper left', framealpha=0.85)

    _add_colorbar(fig, axes.tolist(), cmap, norm, r"Amplitud $A$")

    # Subtítulo con ngap
    fig.suptitle(
        TITLE_ALPHA + rf",  $N_{{\rm gap}} = {NGAP_SIN}$",
        fontsize=13, y=1.01
    )

    out_png = os.path.join(OUT_DIR, "figura_A2_sinusoidal.png")
    
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    
    plt.close()
    print(f"  OK -> {out_png} ")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURA A3 — Gap único + Sinusoidal combinadas 2×3
# ══════════════════════════════════════════════════════════════════════════════

def build_figura_A3():
    print("\n-- FIGURA A3: Gap unico + Sinusoidal 2x3 -----------------------")

    # ── Datos ──────────────────────────────────────────────────────────────
    data_gap = {vf: load_general_r10(vf, ALPHA_CANONICAL) for vf in VFRAGS_A}
    data_sin = {vf: load_sinusoidal_ngap(vf, ALPHA_CANONICAL, NGAP_SIN)
                for vf in VFRAGS_A}

    all_mgaps = sorted(set(e['M_gap'] for vf in VFRAGS_A for e in data_gap[vf]))
    all_amps  = sorted(set(e['amp_val'] for vf in VFRAGS_A for e in data_sin[vf]))
    print(f"  M_gap encontrados: {all_mgaps}")
    print(f"  Amplitudes encontradas: {all_amps}")

    # ── Colormaps ─────────────────────────────────────────────────────────
    cmap_gap = matplotlib.colormaps.get_cmap('inferno')
    norm_gap = LogNorm(vmin=max(min(all_mgaps), 1e-3), vmax=max(all_mgaps))
    def cfn_gap(item): return cmap_gap(norm_gap(item['M_gap']))

    cmap_sin = matplotlib.colormaps.get_cmap('viridis')
    norm_sin = LogNorm(vmin=max(min(all_amps), 0.1), vmax=max(all_amps))
    def cfn_sin(item): return cmap_sin(norm_sin(item['amp_val']))

    # ── Figura 2×3 ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.5),
                             sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.02, hspace=0.06, left=0.07, right=0.85)

    # ── Fila 0: Gap único ─────────────────────────────────────────────────
    for col, vf in enumerate(VFRAGS_A):
        ax = axes[0, col]
        _style_ax(ax)
        entries = data_gap[vf]
        if entries:
            tc = _draw_curves(ax, sorted(entries, key=lambda x: x['M_gap']),
                              cfn_gap)
            _add_snowline(ax, tc)
        else:
            ax.text(0.5, 0.5, 'Sin datos', transform=ax.transAxes,
                    ha='center', va='center', color='gray')
        ax.set_title(rf"$v_{{\rm frag}} = {vf}$ m s$^{{-1}}$",
                     fontsize=17, pad=8)
        if col == 0:
            ax.set_ylabel(YLABEL, fontsize=16)

    # ── Fila 1: Sinusoidal ────────────────────────────────────────────────
    for col, vf in enumerate(VFRAGS_A):
        ax = axes[1, col]
        _style_ax(ax)
        entries = data_sin[vf]
        if entries:
            tc = _draw_curves(ax, sorted(entries, key=lambda x: x['amp_val']),
                              cfn_sin)
            _add_snowline(ax, tc)
        else:
            ax.text(0.5, 0.5, 'Sin datos', transform=ax.transAxes,
                    ha='center', va='center', color='gray')
        ax.set_xlabel(XLABEL, fontsize=16)
        if col == 0:
            ax.set_ylabel(YLABEL, fontsize=16)

    # ── Leyenda común (en panel superior izquierdo) ───────────────────────
    h = [Line2D([0],[0], color='saddlebrown', lw=1.6, ls='--', alpha=0.6),
         Line2D([0],[0], color='black', lw=1.8, ls=':', alpha=0.85),
         Line2D([0],[0], color='royalblue', lw=2.0, ls='--', alpha=0.75)]
    l = [r"$M_{\rm iso}$", r"$0.1\,M_\oplus$", "Snowline"]
    axes[0, 0].legend(h, l, fontsize=14, loc='upper left', framealpha=0.85)

    # ── Colorbar Gap (inferno) — mitad superior derecha ───────────────────
    cbar_ax_gap = fig.add_axes([0.87, 0.535, 0.02, 0.38])
    sm_gap = ScalarMappable(cmap=cmap_gap, norm=norm_gap); sm_gap.set_array([])
    cbar_gap = fig.colorbar(sm_gap, cax=cbar_ax_gap)
    cbar_gap.set_label(r"$M_{\rm gap}\,[M_{\rm Jup}]$", fontsize=16)
    cbar_gap.ax.tick_params(labelsize=15)

    # ── Colorbar Sinusoidal (viridis) — mitad inferior derecha ────────────
    cbar_ax_sin = fig.add_axes([0.87, 0.09, 0.02, 0.38])
    sm_sin = ScalarMappable(cmap=cmap_sin, norm=norm_sin); sm_sin.set_array([])
    cbar_sin = fig.colorbar(sm_sin, cax=cbar_ax_sin)
    cbar_sin.set_label(r"Amplitud $A$", fontsize=16)
    cbar_sin.ax.tick_params(labelsize=15)

    # ── Título común ──────────────────────────────────────────────────────
    fig.suptitle(TITLE_ALPHA, fontsize=18, y=0.98)

    out_png = os.path.join(OUT_DIR, "figura_A3_gap_sinusoidal.png")
    
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    
    plt.close()
    print(f"  OK -> {out_png} ")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURA B — Regímenes M_gap 2×2
# ══════════════════════════════════════════════════════════════════════════════

def build_figura_B():
    print("\n-- FIGURA B: regimenes M_gap 2x2 --------------------------------")

    # Cargar datos por alpha y filtrar por M_gap
    by_mgap = {mg: [] for mg in MGAPS_B}
    for alpha in ALPHAS_B:
        entries = load_general_r10(VFRAG_CANONICAL, alpha)
        for e in entries:
            mg = e.get('M_gap')
            if mg in by_mgap:
                by_mgap[mg].append(e)

    alphas_present = sorted(ALPHAS_B)
    cmap  = matplotlib.colormaps.get_cmap('inferno')
    log_a = np.log10(alphas_present)
    norm  = Normalize(vmin=log_a.min(), vmax=log_a.max())
    def cfn(item): return cmap(norm(np.log10(item['alpha'])))

    fig, axes = plt.subplots(2, 2, figsize=(11, 9),
                             sharex=True, sharey=True)
    axes_flat = axes.flatten()

    for i, mg in enumerate(MGAPS_B):
        ax      = axes_flat[i]
        entries = sorted(by_mgap[mg], key=lambda x: x['alpha'])
        _style_ax(ax)

        # Label del panel: solo M_gap
        ax.text(0.05, 0.97,
                rf"$M_{{\rm gap}} = {mg}\,M_{{\rm Jup}}$",
                transform=ax.transAxes, ha='left', va='top', fontsize=14,
                bbox=dict(facecolor='white', edgecolor='lightgray',
                          alpha=0.9, boxstyle='round,pad=0.3'))

        if entries:
            min_alpha = min(e['alpha'] for e in entries)
            for e in entries:
                color = cfn(e)
                t_myr  = e['times_yr'] / 1e6
                mass_e = e['mass_e']
                mask   = np.isfinite(t_myr) & np.isfinite(mass_e) & (mass_e > 0)
                if not mask.any(): continue
                tv = t_myr[mask]; mv = mass_e[mask]
                if tv[-1] < 10.0:
                    tv = np.append(tv, 10.0); mv = np.append(mv, mv[-1])
                if e['alpha'] == min_alpha and 'm_iso_e' in e:
                    ax.plot(t_myr[mask], e['m_iso_e'][mask],
                            color='saddlebrown', ls='--', alpha=0.5, lw=1.6, zorder=1)
                ax.plot(tv, mv, color=color, alpha=0.88, lw=3.0, zorder=3)

            tc = [e['t_cross_1au']/1e6 for e in entries if e.get('t_cross_1au')]
            _add_snowline(ax, tc)
        else:
            ax.text(0.5, 0.5, 'Sin datos', transform=ax.transAxes,
                    ha='center', va='center', color='gray')

        row, col = i // 2, i % 2
        if col == 0: ax.set_ylabel(YLABEL, fontsize=16)
        if row == 1: ax.set_xlabel(XLABEL, fontsize=16)

    # Leyenda M_iso, Umbral, Snowline en primer panel
    h = [Line2D([0],[0], color='saddlebrown', lw=1.6, ls='--', alpha=0.95),
         Line2D([0],[0], color='black', lw=1.8, ls=':', alpha=0.85),
         Line2D([0],[0], color='royalblue', lw=2.0, ls='--', alpha=0.75)]
    l = [r"$M_{\rm iso}$", r"$0.1\,M_\oplus$", "Snowline"]
    axes_flat[0].legend(h, l, fontsize=13, loc='lower left', framealpha=0.85)

    # Colorbar alpha — labels formateados correctamente
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    fig.subplots_adjust(wspace=0.03, hspace=0.03, right=0.86, top=0.93, bottom=0.09)
    cbar_ax = fig.add_axes([0.88, 0.09, 0.03, 0.84])
    cbar = fig.colorbar(sm, cax=cbar_ax)
    cbar.set_label(r"$\alpha$", fontsize=25, labelpad=15)
    cbar_ticks = np.log10(alphas_present)
    cbar.set_ticks(cbar_ticks)
    def _fmt_alpha(a):
        exp = int(np.floor(np.log10(a) + 1e-9))
        mantissa = a / 10**exp
        if abs(mantissa - 1.0) < 0.05:
            return rf"$10^{{{exp}}}$"
        else:
            return rf"${mantissa:.0f}\times10^{{{exp}}}$"
    cbar.set_ticklabels([_fmt_alpha(a) for a in alphas_present], fontsize=15)

    fig.suptitle(
        rf"$v_{{\rm frag}} = {VFRAG_CANONICAL}$ m s$^{{-1}}$,  "
        rf"$r_{{\rm gap}} = {int(RGAP_FIXED)}$ AU,  "
        r"$M_{\mathrm{emb},0} = 10^{-4}\,M_\oplus$",
        fontsize=14, y=0.975
    )

    out_png = os.path.join(OUT_DIR, "figura_B_regimenes_mgap.png")
    
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    
    plt.close()
    print(f"  OK -> {out_png} ")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURA C — Facet Scatter combinado (General + Sinusoidal) 2×4
# ══════════════════════════════════════════════════════════════════════════════

ALPHAS_FACET = [0.0001, 0.0003, 0.0005, 0.001]

def build_figura_C_facet(x_scale='log', x_max=1.5, suffix=''):
    """Facet scatter 2×4: fila 1 = general, fila 2 = sinusoidal.
    Columnas = alpha. Solo rangos Earth analogs y Waterworlds."""
    import pandas as pd
    print(f"\n-- FIGURA C: Facet Scatter combinado 2x4 {suffix} -----------------------")

    csv_path = os.path.join('data', 'summary_master_todos_casos.csv')
    if not os.path.exists(csv_path):
        print(f"  ERROR: No se encontró {csv_path}")
        return
    df_master = pd.read_csv(csv_path, encoding='utf-8')
    df_master = df_master[(df_master['Is_Valid'] == True) & (np.isclose(df_master['M_emb0'], 1e-3))].copy()
    df_master['f_water'] = df_master['frac_h2o_percent'] / 100.0
    df_master['f_water_plot'] = df_master['f_water'].clip(lower=1e-5)

    scenarios = ['general',"rounded", 'sinusoidal']
    row_labels = ['Gap único', 'Sinusoidal']

    # ── Marcadores por v_frag ─────────────────────────────────────────────
    vfrag_style = {
        1:  {'color': "#3b0f70", 'marker': 'o', 'size': 110},
        3:  {'color': "#c92d59", 'marker': 's', 'size': 120},
        10: {'color': "#f98e09", 'marker': '^', 'size': 130},
    }

    if x_scale == 'log':
        alphas = [0.0001]
        fig_width = 11
        fig_height = 7
    else:
        alphas = ALPHAS_FACET
        fig_width = 20
        fig_height = 9

    # ── Figura 2×ncols ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, len(alphas), figsize=(fig_width, fig_height),
                             sharex=True, sharey=True, squeeze=False)
    fig.subplots_adjust(wspace=0.04, hspace=0.08, left=0.12 if len(alphas)==1 else 0.06, right=0.94 if len(alphas)>1 else 0.88,
                        top=0.82, bottom=0.12)

    np.random.seed(42)

    y_top = max(10, df_master['M_final'].max() * 1.5)
    if np.isnan(y_top) or y_top < 1: y_top = 10.0

    for row in range(2):
        if row == 0:
            df_scen = df_master[df_master['scenario'].str.lower().isin(['general', 'rounded'])]
        else:
            df_scen = df_master[df_master['scenario'].str.lower() == 'sinusoidal']
        for col, alpha in enumerate(alphas):
            ax = axes[row, col]
            sub = df_scen[df_scen['alpha'] == alpha]

            if x_scale == 'log':
                ax.set_xscale('log')
                ax.set_xlim(5e-6, x_max)
            else:
                ax.set_xscale('linear')
                ax.set_xlim(-0.02, 0.52)
            ax.set_yscale('log')
            ax.set_ylim(1e-3, y_top)
            ax.tick_params(which='major', direction='in', top=True, right=True,
                           length=6, width=1.0)
            ax.tick_params(which='minor', direction='in', top=True, right=True,
                           length=3, width=0.7)
            ax.grid(True, which='major', color='lightgray', ls='--', alpha=0.3)

            # Earth analogs: f_H2O 0.01%–0.1% (1e-4 a 1e-3), M 0.1–1 M⊕
            ax.add_patch(plt.Rectangle(
                (1e-4, 0.1), 1e-3 - 1e-4, 1.0 - 0.1,
                facecolor='mediumseagreen', alpha=0.3, zorder=1))
            
            # Super-Earths: f_H2O 0.01%–0.1%, M 1–10 M⊕
            ax.add_patch(plt.Rectangle(
                (1e-4, 1.0), 1e-3 - 1e-4, 10.0 - 1.0,
                facecolor='darkgreen', alpha=0.3, zorder=1))

            # Transición: f_H2O 0.1%–10% (1e-3 a 0.1), M >= 0.1 M⊕
            trans_width = min(0.1, x_max) - 1e-3
            if trans_width > 0:
                ax.add_patch(plt.Rectangle(
                    (1e-3, 0.1), trans_width, y_top - 0.1,
                    facecolor='turquoise', alpha=0.3, zorder=1))

            # Waterworlds: f_H2O ≥ 10%, M ≥ 0.1 M⊕
            ww_width = x_max - 0.1
            if ww_width > 0:
                ax.add_patch(plt.Rectangle(
                    (0.1, 0.1), ww_width, y_top - 0.1,
                    facecolor='dodgerblue', alpha=0.3, zorder=1))

            # Umbral 0.1 M⊕
            ax.axhline(0.1, color='dimgray', ls='--', lw=1.5, alpha=0.8, zorder=2)

            # ── Scatter por v_frag ────────────────────────────────────────
            for vf in sorted(sub['v_frag'].dropna().unique()):
                vf_int = int(vf)
                s = vfrag_style.get(vf_int, {'color': 'gray', 'marker': 'o', 'size': 70})
                sub_vf = sub[sub['v_frag'] == vf]
                jx = np.random.uniform(0.97, 1.03, len(sub_vf))
                jy = np.random.uniform(0.97, 1.03, len(sub_vf))

                label = ''
                if row == 0 and col == 0:
                    label = rf"$v_{{\rm frag}}={vf_int}$ m/s"

                ax.scatter(
                    sub_vf['f_water_plot'].values * jx,
                    sub_vf['M_final'].values * jy,
                    color=s['color'], marker=s['marker'], s=s['size'],
                    alpha=0.9, edgecolors='black', linewidth=0.1,
                    zorder=3, label=label
                )

            # ── Labels de panel ───────────────────────────────────────────
            if row == 0:
                ax.set_title(rf"$\alpha = {alpha:g}$", fontsize=20, pad=12)
            
            # Row labels a la derecha de la última columna
            if col == len(alphas) - 1:
                row_label = "Gap único" if row == 0 else "Sinusoidal"
                ax.text(1.03, 0.5, row_label, transform=ax.transAxes,
                        ha='left', va='center', rotation=-90, fontsize=20)

            # ── Ejes ─────────────────────────────────────────────────────
            ax.set_yticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            ax.set_yticklabels(['-3', '-2', '-1', '0', '1'], fontsize=14)

            if row == 1:
                ax.set_xlabel(r"$f_{\rm H_2O}$", fontsize=24)
                if x_scale == 'log':
                    if x_max <= 0.15:
                        ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1])
                        ax.set_xticklabels(["0%", "0.01%", "0.1%", "1%", "10%"], fontsize=16.5)
                    else:
                        ax.set_xticks([1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 5e-1])
                        ax.set_xticklabels(["0%", "0.01%", "0.1%", "1%", "10%", "50%"], fontsize=16.5)
                else:
                    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
                    ax.set_xticklabels(["0%", "10%", "20%", "30%", "40%", "50%"], fontsize=16.5)
            if col == 0:
                ax.set_ylabel(r"$\log_{10}(M_{\rm emb}\,[M_\oplus])$", fontsize=18)

    # ── Leyenda global ────────────────────────────────────────────────────
    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc='upper center',
                   bbox_to_anchor=(0.5, 0.99), ncol=3, fontsize=20,
                   frameon=False, markerscale=2.4)


    out_png = os.path.join(OUT_DIR, f"figura_C_facet_combinado{suffix}.png")
    
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    
    plt.close()
    print(f"  OK -> {out_png} ")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 60)
    print("  Figuras A1, A2, A3, B, C — §4.2 Crecimiento Individual")
    print("=" * 60)
    build_figura_A1()
    build_figura_A2()
    build_figura_A3()
    build_figura_B()
    build_figura_C_facet(x_scale='linear', x_max=1.05, suffix='_lineal')
    build_figura_C_facet(x_scale='log', x_max=0.1, suffix='_log_10pct')
    print("\nListo.")

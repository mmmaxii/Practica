# -*- coding: utf-8 -*-
"""
analisis_retrasado.py
---------------------
Extrae la evolución del embrión para corridas con gap retrasado.
Genera gráficos de líneas agrupadas por M_gap y t_gap, y mapas de calor.
Las corridas incompletas (<99 snapshots) se grafican con línea punteada.
"""

import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize

plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.titlesize': 18
})

sys.path.append(os.path.abspath(os.path.dirname(__file__)))
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
import dustpy.constants as c

M0_EARTH = 0.01
EMBRYO_AU = 1.0
BASE_DIR = "data/runs/10Myr delayed"
FIG_DIR = f"data/figures/10Myr_retrasado_M0_{M0_EARTH}"
CACHE_FILE = f"data/runs/10Myr_retrasado_cache_M0_{M0_EARTH}.pkl"

os.makedirs(FIG_DIR, exist_ok=True)

def parse_run_name(run_name):
    # e.g., "run_retrasado_r10.0_m0.01_tgap100000"
    parts = run_name.split('_')
    r_val = float(parts[2][1:])
    m_val = float(parts[3][1:])
    tgap_val = float(parts[4][4:])
    return r_val, m_val, tgap_val

def get_target_runs():
    runs = glob.glob(os.path.join(BASE_DIR, "run_*"))
    target_runs = []
    for rpath in runs:
        if os.path.isdir(rpath):
            completed = os.path.exists(os.path.join(rpath, "data0099.hdf5"))
            target_runs.append((rpath, completed))
    target_runs.sort()
    return target_runs

def extract_data():
    runs_info = get_target_runs()
    data = []
    
    print(f"Iniciando extracción PA3 para {len(runs_info)} simulaciones con gap retrasado...")
    for i, (rpath, completed) in enumerate(runs_info):
        run_name = os.path.basename(rpath)
        try:
            r_val, m_val, tgap_val = parse_run_name(run_name)
        except Exception:
            continue
        
        print(f"[{i+1}/{len(runs_info)}] Procesando {run_name} (Completado: {completed})...")
        
        try:
            pa3 = PebbleAccretionModule3.from_datadir(rpath, M_star=1.0)
        except Exception as e:
            print(f"  Error cargando HDF5 en {run_name}: {e}")
            continue
            
        res = pa3.run_growth([EMBRYO_AU], M0_g=M0_EARTH * c.M_earth)
        hist = res[EMBRYO_AU]
        
        if len(hist) == 0:
            print(f"  Sin evolución registrada para {run_name}")
            continue
            
        times_yr = hist[:, 0] / c.year
        mass_e = hist[:, 1] / c.M_earth
        m_h2o = hist[:, 2]
        m_sil = hist[:, 3]
        rsnow_au = hist[:, 4]
        
        total_final = m_h2o[-1] + m_sil[-1]
        frac_h2o_final = 100.0 * (m_h2o[-1] / total_final) if total_final > 0 else 0.0
        
        cross_idx = np.where(rsnow_au < 1.0)[0]
        t_cross_1au = times_yr[cross_idx[0]] if len(cross_idx) > 0 else None
        
        data.append({
            'run_name': run_name,
            'r_gap': r_val,
            'M_gap': m_val,
            't_gap': tgap_val,
            'times_yr': times_yr,
            'mass_e': mass_e,
            'frac_h2o_final': frac_h2o_final,
            't_cross_1au': t_cross_1au,
            'completed': completed
        })
        
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(data, f)
        
    return data

def plot_lines_grouped_by_mgap(data):
    # En este set, r_gap es fijo (10.0). Haremos paneles por M_gap.
    groups = {}
    for item in data:
        mg = item['M_gap']
        if mg not in groups: groups[mg] = []
        groups[mg].append(item)
        
    m_gaps = sorted(list(groups.keys()))
    n_mg = len(m_gaps)
    
    # Grid de subplots
    cols = 3
    rows = (n_mg + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(16, 5 * rows))
    fig.suptitle("Evolución de Masa - Gap Retrasado (R=10 AU)", fontsize=20)
    
    if not isinstance(axes, np.ndarray): axes = np.array([axes])
    axes = axes.flatten()
    
    # Escala de colores para t_gap
    cmap = plt.get_cmap('plasma')
    all_tgaps = sorted(list(set([item['t_gap'] for item in data])))
    if len(all_tgaps) == 0: return
    norm = LogNorm(vmin=min(all_tgaps), vmax=max(all_tgaps))
    
    for i, mg in enumerate(m_gaps):
        ax = axes[i]
        ax.set_title(f"Profundidad: {mg} $M_{{\mathrm{{Jup}}}}$")
        runs_mg = sorted(groups[mg], key=lambda x: x['t_gap'])
        t_cross_avg = []
        
        for item in runs_mg:
            color = cmap(norm(item['t_gap']))
            ls = '-'
            
            t_myr_log = np.log10(item['times_yr'] / 1e6)
            mass_e_log = np.log10(item['mass_e'])
            
            mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
            ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=2.5, ls=ls)
            
            if item['t_cross_1au']: 
                t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
                
        ax.set_xlim([-3, 1])
        
        min_m_local = np.log10(M0_EARTH)
        max_m_local = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in runs_mg])
        padding = (max_m_local - min_m_local) * 0.1 if max_m_local > min_m_local else 0.2
        ax.set_ylim([min_m_local - padding/2, max_m_local + padding])
        
        ax.grid(True, alpha=0.3)
        
        if t_cross_avg:
            t_med = np.median(t_cross_avg)
            label = "Snowline a 1au" if i == 0 else None
            ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0, label=label)
            
        row = i // cols
        col = i % cols
        if col == 0:
            ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
            
        last_in_col = max([idx for idx in range(n_mg) if idx % cols == col])
        if i >= n_mg - cols:
            ax.set_xlabel("Log(t [Myr])")
            
        if i == 0:
            from matplotlib.lines import Line2D
            custom_lines = [Line2D([0], [0], color=cmap(norm(tg)), lw=2.5) for tg in all_tgaps]
            labels = [f"t_gap = {tg/1e6:.1f} Myr" for tg in all_tgaps]
            if t_cross_avg:
                custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
                labels.append("Snowline a 1au")
            ax.legend(custom_lines, labels, loc="upper left", fontsize=10)

    for j in range(n_mg, len(axes)):
        axes[j].set_visible(False)

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f"evo_retrasado_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()

def plot_heatmaps(data):
    all_tgaps = sorted(list(set([d['t_gap'] for d in data])))
    m_gaps = sorted(list(set([d['M_gap'] for d in data])))
    
    Z_mass = np.full((len(m_gaps), len(all_tgaps)), np.nan)
    Z_h2o = np.full((len(m_gaps), len(all_tgaps)), np.nan)
    
    for item in data:
        i = m_gaps.index(item['M_gap'])
        j = all_tgaps.index(item['t_gap'])
        Z_mass[i, j] = item['mass_e'][-1]
        Z_h2o[i, j] = item['frac_h2o_final']
        
    x_idx = np.arange(len(all_tgaps))
    y_idx = np.arange(len(m_gaps))
    
    x_labels = [f"{tg/1e6:.1f}" for tg in all_tgaps]
    
    # 1. Heatmap Masa
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
    mesh1 = ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', 
                      norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
    
    for i in range(len(m_gaps)):
        for j in range(len(all_tgaps)):
            val = Z_mass[i, j]
            if not np.isnan(val):
                text_color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                ax.text(j, i, rf"{val:.2f} $M_\oplus$", ha="center", va="center", color=text_color, fontsize=11)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(x_labels)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(m_gaps)
    
    ax.set_xlabel("Tiempo de formación del Gap ($t_{\mathrm{gap}}$ [Myr])")
    ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
    ax.set_title("Masa Final del Embrión - Gap Retrasado (R=10 AU)")
    plt.savefig(os.path.join(FIG_DIR, f"heatmap_masa_retrasado_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()
    
    # 2. Heatmap H2O
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
    if vmax_h2o == 0: vmax_h2o = 1.0 
    
    mesh2 = ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', 
                      vmin=0, vmax=vmax_h2o)
    
    for i in range(len(m_gaps)):
        for j in range(len(all_tgaps)):
            val = Z_h2o[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=text_color, fontsize=11)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(x_labels)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(m_gaps)
    
    ax.set_xlabel("Tiempo de formación del Gap ($t_{\mathrm{gap}}$ [Myr])")
    ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
    ax.set_title("Fracción de Agua Final - Gap Retrasado (R=10 AU)")
    plt.savefig(os.path.join(FIG_DIR, f"heatmap_h2o_retrasado_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()

def main():
    if os.path.exists(CACHE_FILE):
        print(f"Cargando datos desde caché: {CACHE_FILE}")
        print("Si deseas procesar de nuevo, borra el archivo de caché.")
        with open(CACHE_FILE, 'rb') as f:
            data = pickle.load(f)
    else:
        data = extract_data()
        
    print("\nGenerando gráficos...")
    plot_lines_grouped_by_mgap(data)
    plot_heatmaps(data)
        
    print(f"\n¡Todos los gráficos se han guardado en {FIG_DIR}!")

if __name__ == "__main__":
    main()

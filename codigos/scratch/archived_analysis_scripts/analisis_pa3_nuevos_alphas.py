# -*- coding: utf-8 -*-
"""
analisis_pa3_nuevos_alphas.py
---------------------
Extrae la evolución del embrión para los alphas 0.0005, 0.003 y 0.005, y genera
gráficos de gradiente de evolución temporal y mapas de calor. Las runs que no 
alcanzaron las 99 snapshots se grafican con línea punteada.
"""

import os
import sys
import glob
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
import matplotlib.ticker as ticker

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

M0_EARTH = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0001
EMBRYO_AU = 1.0
BASE_DIR = "data/runs/10Myr"
CACHE_FILE = f"data/runs/10Myr_pa3_cache_nuevos_alphas_M0_{M0_EARTH}.pkl"

TARGET_ALPHAS = [0.0005, 0.003, 0.005]


def parse_run_name(run_name):
    parts = run_name.split('_')
    r_val = float(parts[1][1:])
    m_val = float(parts[2][1:])
    a_val = float(parts[3][1:])
    return r_val, m_val, a_val

def get_target_runs():
    runs = glob.glob(os.path.join(BASE_DIR, "run_*"))
    target_runs = []
    for rpath in runs:
        run_name = os.path.basename(rpath)
        r_gap, M_gap, alpha = parse_run_name(run_name)
        if alpha in TARGET_ALPHAS:
            completed = os.path.exists(os.path.join(rpath, "data0099.hdf5"))
            target_runs.append((rpath, completed))
    target_runs.sort()
    return target_runs

def extract_data():
    runs_info = get_target_runs()
    data = {}
    
    print(f"Iniciando extracción PA3 para {len(runs_info)} simulaciones...")
    for i, (rpath, completed) in enumerate(runs_info):
        run_name = os.path.basename(rpath)
        r_gap, M_gap, alpha = parse_run_name(run_name)
        
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
        
        if alpha not in data:
            data[alpha] = []
            
        data[alpha].append({
            'run_name': run_name,
            'r_gap': r_gap,
            'M_gap': M_gap,
            'times_yr': times_yr,
            'mass_e': mass_e,
            'frac_h2o_final': frac_h2o_final,
            't_cross_1au': t_cross_1au,
            'completed': completed
        })
        
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(data, f)
        
    return data

def plot_lines_grouped_by_rgap(data_alpha, alpha_val, fig_dir):
    groups = {}
    for item in data_alpha:
        rg = item['r_gap']
        if rg not in groups: groups[rg] = []
        groups[rg].append(item)
        
    r_gaps = sorted(list(groups.keys()))
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 11))
    fig.suptitle(rf"$\alpha = {alpha_val}$", fontsize=18)
    
    axes = axes.flatten()
    
    cmap = plt.get_cmap('viridis')
    all_mgaps = sorted(list(set([item['M_gap'] for item in data_alpha])))
    if len(all_mgaps) == 0: return
    norm = LogNorm(vmin=min(all_mgaps), vmax=max(all_mgaps))
    
    for i, (ax, rg) in enumerate(zip(axes, r_gaps)):
        ax.set_title(f"Gap a {rg} AU")
        
        runs_rg = sorted(groups[rg], key=lambda x: x['M_gap'])
        t_cross_avg = []
        
        for item in runs_rg:
            color = cmap(norm(item['M_gap']))
            ls = '-' if item['completed'] else ':'
            
            t_myr_log = np.log10(item['times_yr'] / 1e6)
            mass_e_log = np.log10(item['mass_e'])
            
            mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
            
            ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=2.5, ls=ls)
            if item['t_cross_1au']: t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
            
        ax.set_xlim([-3, 1])
        
        min_m_local = np.nanmin([np.nanmin(np.log10(item['mass_e'])) for item in runs_rg])
        max_m_local = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in runs_rg])
        
        padding = (max_m_local - min_m_local) * 0.1
        if padding == 0: padding = 0.2
        ax.set_ylim([min_m_local - padding, max_m_local + padding])
            
        ax.grid(True, alpha=0.3)
        
        if t_cross_avg:
            t_med = np.median(t_cross_avg)
            label = "Snowline a 1au" if i == 0 else None
            ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0, label=label)
            
        if i % 2 == 0:
            ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
        if i >= 2:
            ax.set_xlabel("Log(t [Myr])")
            
        if i == 0:
            from matplotlib.lines import Line2D
            custom_lines = [Line2D([0], [0], color=cmap(norm(mg)), lw=2.5) for mg in all_mgaps]
            labels = [f"{mg} $M_{{\mathrm{{Jup}}}}$" for mg in all_mgaps]
            custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls=':'))
            labels.append("Incompleta (< 99 snaps)")
            if t_cross_avg:
                custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
                labels.append("Snowline a 1au")
            ax.legend(custom_lines, labels, loc="upper left", fontsize=11)

    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(os.path.join(fig_dir, f"evo_rgap_a{alpha_val}_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()

def plot_lines_grouped_by_mgap(data_alpha, alpha_val, fig_dir):
    groups = {}
    for item in data_alpha:
        mg = item['M_gap']
        if mg not in groups: groups[mg] = []
        groups[mg].append(item)
        
    m_gaps = sorted(list(groups.keys()))
    n_mg = len(m_gaps)
    
    cols = 3
    rows = (n_mg + cols - 1) // cols
    
    fig, axes = plt.subplots(rows, cols, figsize=(16, 5 * rows))
    fig.suptitle(rf"$\alpha = {alpha_val}$", fontsize=18)
    
    if not isinstance(axes, np.ndarray): axes = np.array([axes])
    axes = axes.flatten()
    
    cmap = plt.get_cmap('plasma')
    all_rgaps = sorted(list(set([item['r_gap'] for item in data_alpha])))
    if len(all_rgaps) == 0: return
    norm = Normalize(vmin=min(all_rgaps), vmax=max(all_rgaps))
    
    for i, (ax, mg) in enumerate(zip(axes, m_gaps)):
        ax.set_title(f"Profundidad: {mg} $M_{{\mathrm{{Jup}}}}$")
        runs_mg = sorted(groups[mg], key=lambda x: x['r_gap'])
        t_cross_avg = []
        
        for k, item in enumerate(runs_mg):
            color = cmap(norm(item['r_gap']))
            ls = '-' if item['completed'] else ':'
            
            t_myr_log = np.log10(item['times_yr'] / 1e6)
            mass_e_log = np.log10(item['mass_e'])
            
            mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
            
            if mg == 0.01:
                current_lw = 2.0 + k * 2.0
                z_ord = 10 - k
            else:
                current_lw = 2.5
                z_ord = 2
            
            ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=current_lw, zorder=z_ord, ls=ls)
            if item['t_cross_1au']: t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
            
        ax.set_xlim([-3, 1])
        
        min_m_local = np.nanmin([np.nanmin(np.log10(item['mass_e'])) for item in runs_mg])
        max_m_local = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in runs_mg])
        
        padding = (max_m_local - min_m_local) * 0.1
        if padding == 0: padding = 0.2
        ax.set_ylim([min_m_local - padding, max_m_local + padding])
            
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
        if i == last_in_col:
            ax.set_xlabel("Log(t [Myr])")
            
        if i == 0:
            from matplotlib.lines import Line2D
            custom_lines = [Line2D([0], [0], color=cmap(norm(rg)), lw=2.5) for rg in all_rgaps]
            labels = [f"{rg} AU" for rg in all_rgaps]
            custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls=':'))
            labels.append("Incompleta (< 99 snaps)")
            if t_cross_avg:
                custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
                labels.append("Snowline a 1au")
            ax.legend(custom_lines, labels, loc="upper left", fontsize=11)
            
    for j in range(n_mg, len(axes)):
        axes[j].set_visible(False)
        
    plt.tight_layout()
    plt.subplots_adjust(top=0.92)
    plt.savefig(os.path.join(fig_dir, f"evo_mgap_a{alpha_val}_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()

def plot_heatmaps(data_alpha, alpha_val, fig_dir):
    r_gaps = sorted(list(set([d['r_gap'] for d in data_alpha])))
    m_gaps = sorted(list(set([d['M_gap'] for d in data_alpha])))
    
    Z_mass = np.full((len(m_gaps), len(r_gaps)), np.nan)
    Z_h2o = np.full((len(m_gaps), len(r_gaps)), np.nan)
    
    for item in data_alpha:
        j = r_gaps.index(item['r_gap'])
        i = m_gaps.index(item['M_gap'])
        Z_mass[i, j] = item['mass_e'][-1]
        Z_h2o[i, j] = item['frac_h2o_final']
        
    x_idx = np.arange(len(r_gaps))
    y_idx = np.arange(len(m_gaps))
    
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
    mesh1 = ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', 
                      norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
    
    for i in range(len(m_gaps)):
        for j in range(len(r_gaps)):
            val = Z_mass[i, j]
            if not np.isnan(val):
                text_color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                ax.text(j, i, rf"{val:.3f} $M_\oplus$", ha="center", va="center", color=text_color, fontsize=10)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(r_gaps)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(m_gaps)
    
    ax.set_xlabel("Posición del Gap [AU]")
    ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
    ax.set_title(rf"Masa Final del Embrión ($\alpha = {alpha_val}$)")
    plt.savefig(os.path.join(fig_dir, f"heatmap_masa_a{alpha_val}_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(9, 7))
    vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
    if vmax_h2o == 0: vmax_h2o = 1.0 
    
    mesh2 = ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', 
                      vmin=0, vmax=vmax_h2o)
    
    for i in range(len(m_gaps)):
        for j in range(len(r_gaps)):
            val = Z_h2o[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=text_color, fontsize=11)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(r_gaps)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(m_gaps)
    
    ax.set_xlabel("Posición del Gap [AU]")
    ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
    ax.set_title(rf"Fracción de Agua Final ($\alpha = {alpha_val}$)")
    plt.savefig(os.path.join(fig_dir, f"heatmap_h2o_a{alpha_val}_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()

def main():
    if os.path.exists(CACHE_FILE):
        print(f"Cargando datos desde caché: {CACHE_FILE}")
        print("Si deseas procesar de nuevo, borra el archivo de caché.")
        with open(CACHE_FILE, 'rb') as f:
            data = pickle.load(f)
    else:
        data = extract_data()
        
    print("\nGenerando gráficos para cada alpha...")
    for alpha_val, data_alpha in data.items():
        print(f" -> Generando plots para alpha = {alpha_val} ({len(data_alpha)} corridas)")
        fig_dir = f"data/figures/M_{M0_EARTH}/alpha_{alpha_val}"
        os.makedirs(fig_dir, exist_ok=True)
        plot_lines_grouped_by_rgap(data_alpha, alpha_val, fig_dir)
        plot_lines_grouped_by_mgap(data_alpha, alpha_val, fig_dir)
        plot_heatmaps(data_alpha, alpha_val, fig_dir)
        
    print(f"\n¡Todos los gráficos se han guardado en data/figures/M_{M0_EARTH}!")

if __name__ == "__main__":
    main()

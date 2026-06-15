# -*- coding: utf-8 -*-
"""
analisis_sinusoidal.py
---------------------
Extrae la evolución del embrión para las corridas sinusoidales
(PocosGaps, Intermedio, MuchosGaps) con diferentes amplitudes.
Genera gráficos de líneas agrupadas por tipo de gap y mapas de calor.
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

M0_EARTH = float(sys.argv[1]) if len(sys.argv) > 1 else 0.0001
EMBRYO_AU = 1.0
BASE_DIR = "data/runs/10Myr Sinusoidal"
FIG_DIR = f"data/figures/M_{M0_EARTH}/sinusoidal"
CACHE_FILE = f"data/runs/10Myr_Sinusoidal_cache_M0_{M0_EARTH}.pkl"

os.makedirs(FIG_DIR, exist_ok=True)

# Mapeo de nombres para visualización
GAP_NAMES = {
    "PocosGaps": "Pocos Gaps (tipo HL Tau)",
    "Intermedio": "Intermedio",
    "MuchosGaps": "Muchos Gaps (tipo AS 209)"
}
# Orden deseado en el eje Y del heatmap o subplots
GAP_ORDER = ["PocosGaps", "Intermedio", "MuchosGaps"]

def parse_run_name(run_name):
    # e.g., "PocosGaps_A0.5"
    parts = run_name.split('_')
    gap_type = parts[0]
    amp_val = float(parts[1][1:]) # Ignora la 'A'
    return gap_type, amp_val

def get_target_runs():
    # Busca solo directorios que contengan '_A' en el BASE_DIR
    runs = glob.glob(os.path.join(BASE_DIR, "*_A*"))
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
    
    print(f"Iniciando extracción PA3 para {len(runs_info)} simulaciones sinusoidales...")
    for i, (rpath, completed) in enumerate(runs_info):
        run_name = os.path.basename(rpath)
        try:
            gap_type, amp_val = parse_run_name(run_name)
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
            'gap_type': gap_type,
            'amp_val': amp_val,
            'times_yr': times_yr,
            'mass_e': mass_e,
            'frac_h2o_final': frac_h2o_final,
            't_cross_1au': t_cross_1au,
            'completed': completed
        })
        
    with open(CACHE_FILE, 'wb') as f:
        pickle.dump(data, f)
        
    return data

def plot_lines_grouped_by_gap_type(data):
    # Agrupar datos
    groups = {g: [] for g in GAP_ORDER}
    for item in data:
        if item['gap_type'] in groups:
            groups[item['gap_type']].append(item)
            
    # Amplitudes únicas para la leyenda
    all_amps = sorted(list(set([item['amp_val'] for item in data])))
    if len(all_amps) == 0: return
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle("Evolución de Masa - Discos Sinusoidales", fontsize=20, y=1.05)
    
    cmap = plt.get_cmap('plasma')
    norm = LogNorm(vmin=min(all_amps), vmax=max(all_amps))
    
    for i, gtype in enumerate(GAP_ORDER):
        ax = axes[i]
        title_text = GAP_NAMES.get(gtype, gtype)
        ax.set_title(title_text, pad=10)
        
        runs_gt = sorted(groups[gtype], key=lambda x: x['amp_val'])
        t_cross_avg = []
        
        for item in runs_gt:
            color = cmap(norm(item['amp_val']))
            ls = '-' if item['completed'] else ':'
            
            t_myr_log = np.log10(item['times_yr'] / 1e6)
            mass_e_log = np.log10(item['mass_e'])
            
            mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
            
            ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=2.5, ls=ls)
            
            if item['t_cross_1au']: 
                t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
                
        ax.set_xlim([-3, 1])
        ax.grid(True, alpha=0.3)
        ax.set_xlabel("Log(t [Myr])")
        
        if t_cross_avg:
            t_med = np.median(t_cross_avg)
            label = "Snowline a 1au" if i == 0 else None
            ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0, label=label)
            
        if i == 0:
            ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
            from matplotlib.lines import Line2D
            custom_lines = [Line2D([0], [0], color=cmap(norm(amp)), lw=2.5) for amp in all_amps]
            labels = [f"A = {amp}" for amp in all_amps]
            custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls=':'))
            labels.append("Incompleta (< 99 snaps)")
            if t_cross_avg:
                custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
                labels.append("Snowline a 1au")
            ax.legend(custom_lines, labels, loc="upper left", fontsize=10)
            
    # Ajustar Y limits globales
    min_m_global = np.nanmin([np.nanmin(np.log10(item['mass_e'])) for item in data])
    max_m_global = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in data])
    padding = (max_m_global - min_m_global) * 0.1
    axes[0].set_ylim([min_m_global - padding, max_m_global + padding])

    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, f"evo_sinusoidal_lineas_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()


def plot_heatmaps(data):
    all_amps = sorted(list(set([d['amp_val'] for d in data])))
    
    Z_mass = np.full((len(GAP_ORDER), len(all_amps)), np.nan)
    Z_h2o = np.full((len(GAP_ORDER), len(all_amps)), np.nan)
    
    for item in data:
        if item['gap_type'] not in GAP_ORDER: continue
        i = GAP_ORDER.index(item['gap_type'])
        j = all_amps.index(item['amp_val'])
        
        Z_mass[i, j] = item['mass_e'][-1]
        Z_h2o[i, j] = item['frac_h2o_final']
        
    x_idx = np.arange(len(all_amps))
    y_idx = np.arange(len(GAP_ORDER))
    
    # Nombres bonitos para Y
    y_labels = [GAP_NAMES.get(g, g) for g in GAP_ORDER]
    
    # 1. Heatmap Masa
    fig, ax = plt.subplots(figsize=(10, 6))
    vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
    mesh1 = ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', 
                      norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
    
    for i in range(len(GAP_ORDER)):
        for j in range(len(all_amps)):
            val = Z_mass[i, j]
            if not np.isnan(val):
                text_color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                ax.text(j, i, rf"{val:.2f} $M_\oplus$", ha="center", va="center", color=text_color, fontsize=11)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(all_amps)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(y_labels)
    
    ax.set_xlabel("Amplitud de la perturbación (A)")
    ax.set_title("Masa Final del Embrión - Discos Sinusoidales")
    plt.savefig(os.path.join(FIG_DIR, f"heatmap_masa_sinusoidal_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
    plt.close()
    
    # 2. Heatmap H2O
    fig, ax = plt.subplots(figsize=(10, 6))
    vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
    if vmax_h2o == 0: vmax_h2o = 1.0 
    
    mesh2 = ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', 
                      vmin=0, vmax=vmax_h2o)
    
    for i in range(len(GAP_ORDER)):
        for j in range(len(all_amps)):
            val = Z_h2o[i, j]
            if not np.isnan(val):
                text_color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=text_color, fontsize=11)

    ax.set_xticks(x_idx)
    ax.set_xticklabels(all_amps)
    ax.set_yticks(y_idx)
    ax.set_yticklabels(y_labels)
    
    ax.set_xlabel("Amplitud de la perturbación (A)")
    ax.set_title("Fracción de Agua Final - Discos Sinusoidales")
    plt.savefig(os.path.join(FIG_DIR, f"heatmap_h2o_sinusoidal_M0_{M0_EARTH}.png"), dpi=200, bbox_inches='tight')
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
    plot_lines_grouped_by_gap_type(data)
    plot_heatmaps(data)
        
    print(f"\n¡Todos los gráficos se han guardado en {FIG_DIR}!")

if __name__ == "__main__":
    main()

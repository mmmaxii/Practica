import os, sys, pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Insertamos codigos al path
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))

from pipeline_analysis.analyzer_base import DelayedAnalyzer

# Parámetros por defecto
ALPHA_PLOT = 0.001
RGAP_PLOT = 10.0
MGAP_PLOT = 0.5
M0_EARTH = 1e-4
V_FRAG = 10

OUT_DIR = 'data/figures/paper'
os.makedirs(OUT_DIR, exist_ok=True)

def generate_delayed_cache():
    runs_path = f"data/runs/vf_{V_FRAG}ms/delayed"
    cache_path = os.path.join(runs_path, f"cache_delayed_a_{ALPHA_PLOT}_M0_{M0_EARTH}.pkl")
    
    # Check if cache is empty
    if not os.path.exists(cache_path) or os.path.getsize(cache_path) < 10:
        print(f"Caché vacía o inexistente en {cache_path}. Regenerando...")
        analyzer = DelayedAnalyzer(runs_path, cache_path, M0_EARTH, alpha_val=ALPHA_PLOT)
        return analyzer.load_or_extract()
        
    with open(cache_path, 'rb') as f:
        data = pickle.load(f)
    return data

def plot_resultados_delayed(data, R_gap=RGAP_PLOT, alpha=ALPHA_PLOT):
    unique_mgaps = sorted(list(set([d.get('M_gap', -1) for d in data if np.isclose(d.get('r_gap', -1), R_gap)])))
    mgaps_to_plot = [m for m in unique_mgaps if m > 0]
    available_markers = ['^', 'o', 's', 'D', 'v', '*', 'p', 'h', '<', '>']
    markers = available_markers[:len(mgaps_to_plot)]
    print(f"Generando gráfico 1x1 con markers para M_gap = {mgaps_to_plot}...")
    
    fig, ax = plt.subplots(figsize=(8, 7))
    
    # Colorbar scale limits
    t_gaps_all = [d.get('t_gap', 0.0) / 1e6 for d in data if np.isclose(d.get('r_gap', -1), R_gap)]
    from matplotlib.colors import LogNorm
    vmin = np.min(t_gaps_all) if len(t_gaps_all) > 0 and np.min(t_gaps_all) > 0 else 0.1
    vmax = np.max(t_gaps_all) if len(t_gaps_all) > 0 else 5
    norm = LogNorm(vmin=vmin, vmax=vmax)
    cmap = plt.cm.viridis
    
    # Grid y Línea límite de supervivencia
    ax.grid(True, which='both', linestyle=':', color='black', alpha=0.5)
    ax.axhline(0.1, color='grey', linestyle='--', linewidth=1.5, zorder=2)
    
    # Guardar handles para la leyenda
    legend_handles = []
    
    for mgap, marker in zip(mgaps_to_plot, markers):
        # Filtrar datos
        filtered_data = [d for d in data if np.isclose(d.get('r_gap', -1), R_gap) and np.isclose(d.get('M_gap', -1), mgap)]
            
        frac_h2o = []
        m_final = []
        t_gaps = []
        
        for d in filtered_data:
            mass_arr = d.get('mass_e', [])
            # Filtrar solo simulaciones que llegaron a 98 snapshots
            if len(mass_arr) < 98:
                continue
                
            mf = mass_arr[-1] if len(mass_arr) > 0 else M0_EARTH
            if mf <= 0:
                mf = 1e-4
                
            frac_h2o.append(d.get('frac_h2o_final', 0.0) / 100.0) 
            m_final.append(mf)
            t_gaps.append(d.get('t_gap', 0.0) / 1e6)
            
        if len(t_gaps) > 0:
            sc = ax.scatter(frac_h2o, m_final, c=t_gaps, cmap=cmap, norm=norm, marker=marker, s=180, edgecolors='black', linewidths=1.0, zorder=3)
            # Para la leyenda, creamos un punto negro con el marker correspondiente
            legend_handles.append(plt.Line2D([0], [0], marker=marker, color='w', markerfacecolor='gray', markeredgecolor='black', markersize=12, label=rf'${mgap}\,M_{{\rm Jup}}$'))

    if legend_handles:
        ax.legend(handles=legend_handles, title_fontsize=18, 
                 fontsize=18, loc='upper left', framealpha=0.4)
    
    # Configuración de ejes
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1e-1)
    ax.set_xlim(0, 0.50)
    
    ax.set_xlabel(r'$f_{\rm H_2O}$', fontsize=20)
    ax.set_ylabel(r'$\log_{10}(M_{\rm emb}\,[M_\oplus])$', fontsize=20)
    
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: f"{int(x*100)}%"))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, pos: f"{int(np.log10(y))}"))
    ax.tick_params(axis='both', which='major', labelsize=18)
    
    # Título principal
    ax.set_title(rf"$\alpha={alpha}$, $R_{{\rm gap}}={R_gap}$ au", fontsize=24, pad=15)
    
    # Ajustar espacio para el colorbar
    plt.tight_layout()
    
    # Colorbar
    cbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax, pad=0.02)
    cbar.set_label(r'$t_{\rm gap}\,[\rm Myr]$', fontsize=20)
    
    cbar_ticks = [0.1, 0.5, 1.0, 2.0, 5.0]
    cbar_ticklabels = ['0.1', '0.5', '1', '2', '5']
    cbar.set_ticks(cbar_ticks)
    cbar.set_ticklabels(cbar_ticklabels)
    
    cbar.ax.tick_params(labelsize=18)
    
    out_file = os.path.join(OUT_DIR, f"resultados_delayed_1x1_markers.png")
    fig.savefig(out_file, dpi=350, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura guardada en {out_file}")

if __name__ == '__main__':
    data = generate_delayed_cache()
    plot_resultados_delayed(data)

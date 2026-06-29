import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "axes.labelsize": 18,
    "font.size": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

def main():
    print("Generando Bubble Chart 3x3...")
    
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    
    # Filtro: M_emb0 = 1e-4
    df = df[np.isclose(df['M_emb0'], 0.001)].copy()
    
    v_frags = [1, 3, 10]
    alphas = [0.0001, 0.0005, 0.001]
    
    plasma = cm.get_cmap('plasma')
    cmap_custom = LinearSegmentedColormap.from_list('trunc_plasma', plasma(np.linspace(0.0, 0.85, 256)))
    
    fig, axes = plt.subplots(3, 3, figsize=(14, 14), sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0.04, hspace=0.04, right=0.825)
    
    for row, v_frag in enumerate(v_frags):
        for col, alpha in enumerate(alphas):
            ax = axes[row, col]
            
            df_plot = df[(df['v_frag'] == v_frag) & np.isclose(df['alpha'], alpha)]
            
            ax.tick_params(direction='in', top=True, right=True, which='both', labelsize=18)
            ax.grid(True, which='major', color='lightgray', ls='--', alpha=0.3)

            # Titles on top row
            if row == 0:
                ax.set_title(rf"$\alpha = {alpha:g}$", fontsize=20, pad=12)
            
            # Row labels on the right side of the last column
            if col == 2:
                ax.text(1.05, 0.5, rf"$v_{{\rm frag}} = {v_frag}$ m/s", 
                        transform=ax.transAxes, rotation=-90, va='center', ha='left', fontsize=20)
            
            if row == 0 and col == 2:
                legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                               r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
                ax.text(0.95, 0.95, legend_text,
                        transform=ax.transAxes, fontsize=20,
                        ha='right', va='top',
                        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
                        
            if len(df_plot) == 0:
                continue
                
            mask_success = df_plot['M_final'] >= 0.1
            mask_fail = df_plot['M_final'] < 0.1
            
            df_succ = df_plot[mask_success]
            df_fail = df_plot[mask_fail]
            
            if len(df_fail) > 0:
                ax.scatter(
                    df_fail['r_gap'], df_fail['M_gap'],
                    s=60, c='#aaaaaa', marker='v', alpha=0.7,
                    edgecolors='black', linewidth=0.5
                )
                
            if len(df_succ) > 0:
                sizes = df_succ['M_final'] * 200 + 20
                frac_h2o = df_succ['frac_h2o_percent'] / 100.0
                
                sc = ax.scatter(
                    df_succ['r_gap'], df_succ['M_gap'],
                    s=sizes, c=frac_h2o,
                    cmap=cmap_custom, vmin=0.0, vmax=0.10,
                    marker='o', alpha=0.85,
                    edgecolors='black', linewidth=0.7
                )
            
            ax.set_yscale('log')
            ax.set_ylim(0.005, 5.0)
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.set_yticks([0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 5.0])
            ax.set_xlim(0, 32)
            ax.set_xticks([0, 5, 10, 15, 20, 25, 30])
            
            if col == 0:
                ax.set_ylabel(r"$M_{\rm gap}$ [$M_{\rm Jup}$]", fontsize=22)
            if row == 2:
                ax.set_xlabel(r"$r_{\rm gap}$ [AU]", fontsize=22)

    # Meta-information (Colorbar) on the right side
    cax = fig.add_axes([0.88, 0.15, 0.03, 0.70])
    sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.10))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
    cbar.set_label("Fracción de Agua Final", fontsize=25)
    cbar.ax.tick_params(labelsize=20)
    cbar.set_ticks([0.0, 0.05, 0.1])
    cbar.set_ticklabels(['0%', '5%', r'$\geq 10\%$'])
    
    out_dir = "data/figures/paper"
    os.makedirs(out_dir, exist_ok=True)
    out_path_png = os.path.join(out_dir, "resultados_bubble.png")
    plt.savefig(out_path_png, dpi=350, bbox_inches='tight')
    plt.close(fig)
    print(f"Completado. Guardado en {out_path_png}")

if __name__ == '__main__':
    main()

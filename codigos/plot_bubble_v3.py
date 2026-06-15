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
    "axes.labelsize": 12,
    "font.size": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

def main():
    print("Generando Bubble Chart 2x3 para v_frag = 3 m/s...")
    
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    
    # Filtros: v_frag = 3, M_emb0 = 0.001 (baseline)
    df = df[(df['v_frag'] == 3) & np.isclose(df['M_emb0'], 0.001)].copy()
    
    # alphas solicitados
    alphas = [0.0001, 0.0003, 0.0005, 0.0007, 0.001]
    
    # Colormap
    custom_colors = ['#d73027', '#0077bb', '#660099'] 
    cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
    
    # Grilla 2x3
    fig, axes = plt.subplots(2, 3, figsize=(14, 9), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    axes_flat = axes.flatten()
    
    for i, alpha in enumerate(alphas):
        ax = axes_flat[i]
        df_alpha = df[np.isclose(df['alpha'], alpha)]
        
        ax.tick_params(direction='in', top=True, right=True, which='both')
        
        if len(df_alpha) == 0:
            ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='center', va='top', fontsize=13, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
            continue
            
        mask_success = df_alpha['M_final'] >= 0.1
        mask_fail = df_alpha['M_final'] < 0.1
        
        df_succ = df_alpha[mask_success]
        df_fail = df_alpha[mask_fail]
        
        # Triángulos (Fracasos)
        if len(df_fail) > 0:
            ax.scatter(
                df_fail['r_gap'], df_fail['M_gap'],
                s=40, c='#d3d3d3', marker='v', alpha=0.7,
                edgecolors='black', linewidth=0.5
            )
            
        # Círculos (Éxitos)
        if len(df_succ) > 0:
            # Tamaño proporcional a M_final
            sizes = df_succ['M_final'] * 400 + 20
            # Fracción de agua de 0 a 1 -> lo pasamos a escala 0 a 0.20 para el cmap
            frac_h2o = df_succ['frac_h2o_percent'] / 100.0
            
            sc = ax.scatter(
                df_succ['r_gap'], df_succ['M_gap'],
                s=sizes, c=frac_h2o,
                cmap=cmap_custom, vmin=0.0, vmax=0.20,
                marker='o', alpha=0.85,
                edgecolors='black', linewidth=0.7
            )
            
        # Título interno
        ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                ha='center', va='top', fontsize=13, 
                bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
        
        # Escala logarítmica para Y
        ax.set_yscale('log')
        ax.set_ylim(0.005, 5.0)
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.set_yticks([0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 5.0])
        ax.set_xlim(min(df['r_gap'])-2, max(df['r_gap'])+2)

        # Labels solo en bordes izquierdos y bajos
        if i % 3 == 0:
            ax.set_ylabel(r"Profundidad del Gap [$M_{\rm Jup}$]")
        if i >= 3:
            ax.set_xlabel("Posición del Gap [AU]")

    # ====================================================
    # 6TO PANEL: METAINFORMACIÓN
    # ====================================================
    ax_meta = axes_flat[5]
    ax_meta.axis('off')  # Ocultar ejes
    
    # 1. Colorbar (Water fraction)
    # Insertamos un inset_axes dentro del ax_meta para la barra
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    cax = inset_axes(ax_meta, width="70%", height="10%", loc='upper center', bbox_to_anchor=(0, -0.15, 1, 1), bbox_transform=ax_meta.transAxes)
    sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
    sm.set_array([])
    cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label("Fracción de Agua Final")
    cbar.set_ticks([0.0, 0.1, 0.2])
    cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
    
    # 2. Leyenda de Éxito / Fracaso
    legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                   r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
    ax_meta.text(0.5, 0.5, legend_text,
            transform=ax_meta.transAxes, fontsize=13,
            ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
            
    # 3. Leyenda de Tamaños de Burbuja
    ax_meta.text(0.5, 0.2, "Tamaño: Masa Final ($M_\oplus$)",
            transform=ax_meta.transAxes, fontsize=11, ha='center', va='center')
    ax_meta.scatter([0.3], [0.1], s=0.1 * 400 + 20, color='gray', transform=ax_meta.transAxes)
    ax_meta.text(0.3, 0.05, "0.1", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')
    ax_meta.scatter([0.5], [0.1], s=0.5 * 400 + 20, color='gray', transform=ax_meta.transAxes)
    ax_meta.text(0.5, 0.05, "0.5", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')
    ax_meta.scatter([0.7], [0.1], s=1.0 * 400 + 20, color='gray', transform=ax_meta.transAxes)
    ax_meta.text(0.7, 0.05, "1.0", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')

    # Ajuste final
    fig.suptitle(r"Formación de Waterworlds a $v_{\rm frag} = 3$ m/s ($M_{\rm emb, 0} = 10^{-3} M_\oplus$)", fontsize=16, y=0.92)
    
    out_path = "data/benchmarks/bubble_mosaic_v3.png"
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Completado. Guardado en {out_path}")

if __name__ == '__main__':
    main()

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import ScalarFormatter

from pipeline_analysis.snapshot_analyzer import SnapshotAnalyzer
from pipeline_analysis.plotters import PlotterBenchmarks

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

DIRECTORIO_SALIDA = "data/benchmarks/"
os.makedirs(DIRECTORIO_SALIDA, exist_ok=True)

def generar_bubble_mosaic(v_frag_target):
    print(f"\n[1/4] Generando Bubble Chart Mosaico para v_frag = {v_frag_target} m/s...")
    
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"   [!] Error: No se encontró {csv_path}. Ejecuta compilar_tabla_global.py primero.")
        return
        
    df = pd.read_csv(csv_path)
    
    # Filtros: v_frag seleccionado
    # IMPORTANTE: El usuario ha confirmado que M_emb0 siempre es 1e-4 (0.0001)
    target_m0 = 0.0001
        
    df = df[(df['v_frag'] == v_frag_target) & np.isclose(df['M_emb0'], target_m0)].copy()
    
    # Alphas estrictos para asegurar consistencia
    if v_frag_target in [1, 3]:
        alphas = [0.0001, 0.0003, 0.0005, 0.0007, 0.001]
    else:
        alphas = [0.0001, 0.0005, 0.001, 0.003, 0.005, 0.01]
    
    if len(df) == 0:
        print(f"   [!] No hay datos en la tabla maestra para v_frag = {v_frag_target}")
        return
        
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
            ax.set_title(rf"$\alpha = {alpha}$", fontsize=13)
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
            # Tamaño proporcional a M_final reducido para evitar superposiciones masivas
            factor_escala = 100
            sizes = df_succ['M_final'] * factor_escala + 20
            # Fracción de agua
            frac_h2o = df_succ['frac_h2o_percent'] / 100.0
            
            sc = ax.scatter(
                df_succ['r_gap'], df_succ['M_gap'],
                s=sizes, c=frac_h2o,
                cmap=cmap_custom, vmin=0.0, vmax=0.20,
                marker='o', alpha=0.75,
                edgecolors='black', linewidth=0.5
            )
            
        # Título externo para que no choque con los puntos altos
        ax.set_title(rf"$\alpha = {alpha}$", fontsize=13)
        
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

    # PANEL DE METAINFORMACIÓN DEPENDIENDO DE LA CANTIDAD DE ALPHAS
    if len(alphas) == 5:
        # 6TO PANEL LIBRE: Usamos el diseño con leyenda completa de tamaños
        ax_meta = axes_flat[5]
        ax_meta.axis('off')  
        
        from mpl_toolkits.axes_grid1.inset_locator import inset_axes
        cax = inset_axes(ax_meta, width="70%", height="10%", loc='upper center', bbox_to_anchor=(0, -0.15, 1, 1), bbox_transform=ax_meta.transAxes)
        sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cax, orientation='horizontal')
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        ax_meta.text(0.5, 0.6, legend_text,
                transform=ax_meta.transAxes, fontsize=13,
                ha='center', va='center',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
                
        factor_escala = 100
        ax_meta.text(0.5, 0.3, "Tamaño: Masa Final ($M_\oplus$)",
                transform=ax_meta.transAxes, fontsize=11, ha='center', va='center')
        ax_meta.scatter([0.3], [0.15], s=0.1 * factor_escala + 20, color='gray', transform=ax_meta.transAxes)
        ax_meta.text(0.3, 0.05, "0.1", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')
        ax_meta.scatter([0.5], [0.15], s=1.0 * factor_escala + 20, color='gray', transform=ax_meta.transAxes)
        ax_meta.text(0.5, 0.05, "1.0", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')
        ax_meta.scatter([0.7], [0.15], s=3.0 * factor_escala + 20, color='gray', transform=ax_meta.transAxes)
        ax_meta.text(0.7, 0.05, "3.0", transform=ax_meta.transAxes, fontsize=10, ha='center', va='center')
    else:
        # DISEÑO DE 6 PANELES OCUPADOS: Colorbar a la derecha, leyenda en panel 0
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        axes_flat[5].text(0.05, 0.15, legend_text,
                transform=axes_flat[5].transAxes, fontsize=11,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        
        factor_escala = 100
        axes_flat[5].text(0.65, 0.35, "Tamaño: M. Final ($M_\oplus$)",
                transform=axes_flat[5].transAxes, fontsize=9, ha='center', va='center')
        axes_flat[5].scatter([0.45], [0.2], s=0.1 * factor_escala + 20, color='gray', transform=axes_flat[5].transAxes)
        axes_flat[5].text(0.45, 0.1, "0.1", transform=axes_flat[5].transAxes, fontsize=9, ha='center', va='center')
        axes_flat[5].scatter([0.65], [0.2], s=1.0 * factor_escala + 20, color='gray', transform=axes_flat[5].transAxes)
        axes_flat[5].text(0.65, 0.1, "1.0", transform=axes_flat[5].transAxes, fontsize=9, ha='center', va='center')
        axes_flat[5].scatter([0.85], [0.2], s=3.0 * factor_escala + 20, color='gray', transform=axes_flat[5].transAxes)
        axes_flat[5].text(0.85, 0.1, "3.0", transform=axes_flat[5].transAxes, fontsize=9, ha='center', va='center')

    fig.suptitle(rf"Formación de Waterworlds a $v_{{\rm frag}} = {v_frag_target}$ m/s ($M_{{\rm emb, 0}} = 10^{{-4}} M_\oplus$)", fontsize=16, y=0.92)
    
    out_dir = os.path.join(DIRECTORIO_SALIDA, f"vf_{v_frag_target}ms")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, f"bubble_mosaic_v{v_frag_target}.png")
    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f" -> Guardado: {out_path}")

def generar_espaciotemporales(v_frag_target):
    # Alphas estrictos garantizados
    if v_frag_target in [1, 3]:
        alphas_mosaic = ['0.0001', '0.0003', '0.0005', '0.0007', '0.001']
    else:
        alphas_mosaic = ['0.0001', '0.0005', '0.001', '0.003', '0.005', '0.01']

    data_hov = {}

    print(f"\n[2/4] Extrayendo datos espaciotemporales para v_frag = {v_frag_target} m/s...")
    print(f"   -> Alphas a analizar: {alphas_mosaic}")
    
    for alpha in alphas_mosaic:
        # Ajustar el nombre del directorio y archivo
        suffix = f"_v{float(v_frag_target)}" if v_frag_target in [1, 3] else ""
        
        # Soportamos tanto carpetas 'general' como 'rounded'
        arq_path_gen = f"data/runs/vf_{v_frag_target}ms/general/run_r10.0_m0.3_a{alpha}{suffix}"
        arq_path_rnd = f"data/runs/vf_{v_frag_target}ms/rounded/run_r10.0_m0.3_a{alpha}{suffix}"
        
        if os.path.exists(arq_path_gen):
            arq_path = arq_path_gen
        elif os.path.exists(arq_path_rnd):
            arq_path = arq_path_rnd
        else:
            arq_path = None
        
        if arq_path:
            print(f"   -> Procesando SnapshotAnalyzer: {arq_path}")
            analyzer = SnapshotAnalyzer(arq_path)
            t_arr, r_arr, eps_mat, flux_arr, amax_arr = analyzer.extract_spatiotemporal_data(subsampling=1)
            data_hov[alpha] = {
                't_arr': t_arr,
                'r_arr': r_arr,
                'eps_mat': eps_mat,
                'flux_arr': flux_arr,
                'amax_arr': amax_arr
            }
        else:
            print(f"   [!] Faltan datos crudos (ni general ni rounded) para alpha={alpha}")

    if data_hov:
        print("\n[3/4] Generando mosaico de Hovmöller...")
        title_str = rf"Evolución Espaciotemporal a $v_{{\rm frag}} = {v_frag_target}$ m/s"
        out_dir = os.path.join(DIRECTORIO_SALIDA, f"vf_{v_frag_target}ms")
        os.makedirs(out_dir, exist_ok=True)
        out_hov = os.path.join(out_dir, f"hovmoller_mosaico_alphas_v{v_frag_target}.png")
        PlotterBenchmarks.plot_hovmoller_mosaic(data_hov, alphas_mosaic, out_hov, title=title_str)
        print(f" -> Guardado: {out_hov}")
        
        print("\n[4/4] Generando mosaico de Pebble Flux + A_max...")
        out_amax = os.path.join(out_dir, f"pebble_flux_amax_mosaico_v{v_frag_target}.png")
        PlotterBenchmarks.plot_flux_amax_mosaic(data_hov, alphas_mosaic, out_amax, title=title_str)
        print(f" -> Guardado: {out_amax}")
    else:
        print("\n[!] No se pudo encontrar ningún dato para generar los mosaicos espaciotemporales.")

def generar_espaciotemporales_smooth(v_frag_target):
    if v_frag_target in [1, 3]:
        alphas_mosaic = ['0.0001', '0.0003', '0.0005', '0.0007', '0.001']
    else:
        alphas_mosaic = ['0.0001', '0.0005', '0.001', '0.003', '0.005', '0.01']

    data_hov = {}
    print(f"\n[X] Extrayendo datos espaciotemporales Smooth para v_frag = {v_frag_target} m/s...")
    
    import glob
    
    for alpha in alphas_mosaic:
        base_dir = f"data/runs/vf_{v_frag_target}ms/smooth"
        pattern = os.path.join(base_dir, f"run_smooth_a{alpha}*")
        matches = glob.glob(pattern)
        
        if matches:
            arq_path = matches[0]
            print(f"   -> Procesando: {arq_path}")
            analyzer = SnapshotAnalyzer(arq_path)
            t_arr, r_arr, eps_mat, flux_arr, amax_arr = analyzer.extract_spatiotemporal_data(subsampling=1)
            data_hov[alpha] = {
                't_arr': t_arr, 'r_arr': r_arr, 'eps_mat': eps_mat,
                'flux_arr': flux_arr, 'amax_arr': amax_arr
            }
            
    if data_hov:
        out_dir = os.path.join(DIRECTORIO_SALIDA, f"vf_{v_frag_target}ms")
        os.makedirs(out_dir, exist_ok=True)
        title_str = rf"Evolución Espaciotemporal Smooth ($v_{{\rm frag}} = {v_frag_target}$ m/s)"
        
        out_hov = os.path.join(out_dir, f"hovmoller_mosaico_smooth_v{v_frag_target}.png")
        PlotterBenchmarks.plot_hovmoller_mosaic(data_hov, alphas_mosaic, out_hov, title=title_str)
        print(f" -> Guardado: {out_hov}")
        
        out_amax = os.path.join(out_dir, f"pebble_flux_amax_smooth_v{v_frag_target}.png")
        PlotterBenchmarks.plot_flux_amax_mosaic(data_hov, alphas_mosaic, out_amax, title=title_str)
        print(f" -> Guardado: {out_amax}")
    else:
        print("\n[!] No se encontraron datos Smooth para generar los mosaicos espaciotemporales.")

def generar_espaciotemporales_sinusoidal(v_frag_target):
    import glob
    
    import glob
    
    base_dir_search = f"data/runs/vf_{v_frag_target}ms/sinusoidal"
    if not os.path.exists(base_dir_search): return
    
    # Extraer alphas únicos que realmente existen en las carpetas
    all_runs = glob.glob(os.path.join(base_dir_search, "run_ngap*"))
    alphas_encontrados = set()
    import re
    for run_path in all_runs:
        m = re.search(r'_a([0-9.]+)', run_path)
        if m:
            alphas_encontrados.add(m.group(1))
            
    if not alphas_encontrados:
        print(f"   [!] No se encontraron carpetas con alphas para Sinusoidal a {v_frag_target} m/s")
        return
        
    alphas_mosaic = sorted(list(alphas_encontrados), key=float)
        
    amps_mosaic = ['0.5', '0.7', '1.0', '2.0', '3.0', '5.0']
    ngaps_list = [3, 5, 10, 15, 20]
    
    print(f"\n[X] Extrayendo datos espaciotemporales Sinusoidal para v_frag = {v_frag_target} m/s...")
    
    for ngap in ngaps_list:
        for alpha in alphas_mosaic:
            data_hov = {}
            for amp in amps_mosaic:
                base_dir = f"data/runs/vf_{v_frag_target}ms/sinusoidal"
                if not os.path.exists(base_dir): continue
                
                pattern = os.path.join(base_dir, f"run_ngap{ngap}_A{amp}_a{alpha}*")
                matches = glob.glob(pattern)
                if matches:
                    arq_path = matches[0]
                    print(f"   -> Procesando: {arq_path}")
                    analyzer = SnapshotAnalyzer(arq_path)
                    t_arr, r_arr, eps_mat, flux_arr, amax_arr = analyzer.extract_spatiotemporal_data(subsampling=1)
                    data_hov[amp] = {
                        't_arr': t_arr, 'r_arr': r_arr, 'eps_mat': eps_mat,
                        'flux_arr': flux_arr, 'amax_arr': amax_arr
                    }
                    
            if data_hov:
                out_dir = os.path.join(DIRECTORIO_SALIDA, f"vf_{v_frag_target}ms", "hov_sin")
                os.makedirs(out_dir, exist_ok=True)
                title_str = rf"Evolución Espaciotemporal Sinusoidal ($v_{{\rm frag}} = {v_frag_target}$ m/s, $\alpha={alpha}$, $N_{{\rm gaps}}={ngap}$)"
                
                out_hov = os.path.join(out_dir, f"hovmoller_mosaico_ngap{ngap}_a{alpha}.png")
                PlotterBenchmarks.plot_hovmoller_mosaic(data_hov, amps_mosaic, out_hov, title=title_str, r_min_val=0.7)
                print(f" -> Guardado: {out_hov}")
                
                out_amax = os.path.join(out_dir, f"pebble_flux_amax_ngap{ngap}_a{alpha}.png")
                PlotterBenchmarks.plot_flux_amax_mosaic(data_hov, amps_mosaic, out_amax, title=title_str)
                print(f" -> Guardado: {out_amax}")
            else:
                pass # print(f"   [!] No se encontraron datos para ngap={ngap}, alpha={alpha}")

def generar_benchmark_smooth():
    print(f"\nGenerando Benchmark Global Smooth...")
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"   [!] Error: No se encontró {csv_path}.")
        return
        
    df = pd.read_csv(csv_path)
    target_m0 = 0.0001
    
    df_smooth = df[(df['scenario'] == 'smooth') & np.isclose(df['M_emb0'], target_m0)].copy()
    
    if len(df_smooth) == 0:
        print(f"   [!] No hay datos en la tabla maestra para Smooth")
        return
        
    out_path = os.path.join(DIRECTORIO_SALIDA, "bubble_mosaic_smooth.png")
    PlotterBenchmarks.plot_smooth_global_benchmark(df_smooth, out_path)
    print(f" -> Guardado: {out_path}")

def generar_bubble_mosaic_sinusoidal(v_frag_target):
    print(f"\nGenerando Bubble Mosaic Sinusoidal para v_frag = {v_frag_target} m/s...")
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path): return
    df = pd.read_csv(csv_path)
    
    # Filtrar por sinusoidal, v_frag y M_emb0 = 0.0001 para evitar superposiciones
    target_m0 = 0.0001
    df_sin = df[(df['scenario'] == 'sinusoidal') & (df['v_frag'] == v_frag_target) & np.isclose(df['M_emb0'], target_m0)]
    if len(df_sin) == 0:
        print(" [!] No hay datos sinusoidales para esta velocidad y M_emb0.")
        return
        
    alphas_mosaic = sorted(df_sin['alpha'].unique())
        
    data_por_alpha = {}
    for alpha in alphas_mosaic:
        df_alpha = df_sin[np.isclose(df_sin['alpha'], alpha, rtol=1e-5)]
        data_por_alpha[alpha] = {
            'amp_val': df_alpha['amp_val'].values,
            'n_gaps': df_alpha['n_gaps'].values,
            'M_final': df_alpha['M_final'].values,
            'frac_h2o': df_alpha['frac_h2o_percent'].values / 100.0
        }
        
    out_dir = os.path.join(DIRECTORIO_SALIDA, f"vf_{v_frag_target}ms")
    os.makedirs(out_dir, exist_ok=True)
    out_file = os.path.join(out_dir, f"bubble_mosaic_sinusoidal_v{v_frag_target}.png")
    
    PlotterBenchmarks.plot_sinusoidal_bubble_mosaic(data_por_alpha, alphas_mosaic, out_file)
    print(f" -> Guardado: {out_file}")

def main():
    print("======================================================")
    print("   GENERADOR UNIFICADO DE MOSAICOS PARA MANUSCRITO    ")
    print("======================================================")
    v_frag_input = input("Ingrese v_frag a analizar/graficar (ej: 1, 3, 10) o 'todos' para generar todo: ").strip().lower()
    
    if v_frag_input == 'todos':
        v_frags = [1, 3, 10]
    else:
        try:
            v_frags = [int(v_frag_input)]
        except:
            print("Entrada inválida. Usando v_frag = 1 por defecto.")
            v_frags = [1]

    print("Opciones:")
    print("[1] Generar todo secuencialmente (General, Smooth, Sinusoidal)")
    print("[2] Solo General/Rounded (por v_frag)")
    print("[3] Solo Smooth (Global Bubble + Hovmoller por v_frag)")
    print("[4] Solo Sinusoidal (Global Bubble + Hovmoller por v_frag)")
    choice = input("Elija una opción [1]: ").strip() or "1"
    
    # Benchmarks globales (se hacen una sola vez)
    if choice in ['1', '3']:
        generar_benchmark_smooth()
        
        
    for v_frag_target in v_frags:
        print(f"\n--- Iniciando generación de Mosaicos Espaciotemporales (v_frag target = {v_frag_target} m/s) ---")
        
        if choice in ['1', '2']:
            generar_bubble_mosaic(v_frag_target)
            generar_espaciotemporales(v_frag_target)
            
        if choice in ['1', '3']:
            generar_espaciotemporales_smooth(v_frag_target)
            
        if choice in ['1', '4']:
            generar_bubble_mosaic_sinusoidal(v_frag_target)
            # generar_espaciotemporales_sinusoidal(v_frag_target)

    print(f"\n[ÉXITO] Todos los mosaicos generados en {DIRECTORIO_SALIDA}")

if __name__ == '__main__':
    main()

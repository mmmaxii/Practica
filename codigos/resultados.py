import os
import glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import FuncFormatter
from pipeline_analysis.snapshot_analyzer import SnapshotAnalyzer
from pipeline_methods.oka_interpolation import r_snow_time_cgs
import dustpy.constants as c

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "axes.labelsize": 16,
    "font.size": 14,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
})

DIRECTORIO_SALIDA = "data/figures/resultados"
os.makedirs(DIRECTORIO_SALIDA, exist_ok=True)

def find_run_path(scenario, v_frag, alpha, m_gap=None, n_gap=None, amp=None):
    base_dir = f"data/runs/vf_{v_frag}ms"
    if scenario == "smooth":
        pattern = os.path.join(base_dir, "smooth", f"run_smooth_a{alpha}*")
    elif scenario == "sinusoidal":
        pattern = os.path.join(base_dir, "sinusoidal", f"run_ngap{n_gap}_A{amp}_a{alpha}*")
        matches = glob.glob(pattern)
        if matches:
            return matches[0]
        return None
    else:
        pattern_gen = os.path.join(base_dir, "general", f"run_r10.0_m{m_gap}_a{alpha}*")
        pattern_rnd = os.path.join(base_dir, "rounded", f"run_r10.0_m{m_gap}_a{alpha}*")
        
        matches = glob.glob(pattern_gen)
        if not matches:
            matches = glob.glob(pattern_rnd)
        if matches:
            return matches[0]
        return None
        
    matches = glob.glob(pattern)
    if matches:
        return matches[0]
    return None

def extract_data(path):
    if not path or not os.path.exists(path):
        print(f"   [!] Path no encontrado: {path}")
        return None
    print(f"   -> Extrayendo: {path}")
    analyzer = SnapshotAnalyzer(path)
    t_arr, r_arr, eps_mat, _, _ = analyzer.extract_spatiotemporal_data(subsampling=1)
    return {'t_arr': t_arr, 'r_arr': r_arr, 'eps_mat': eps_mat}

def plot_hovmoller_grid(data_list, titles, out_name, rows, cols, figsize=(12, 8), suptitle=None, col_titles=None):
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.ticker import FuncFormatter
    import numpy as np
    
    # Check what variables are global
    global r_snow_time_cgs
    
    fig, axes = plt.subplots(rows, cols, figsize=figsize, sharex=True, sharey=True)
    if rows == 1 and cols == 1:
        axes = np.array([[axes]])
    elif rows == 1:
        axes = axes[np.newaxis, :]
    elif cols == 1:
        axes = axes[:, np.newaxis]
        
    plt.subplots_adjust(wspace=0.0, hspace=0.0)
    
    eps_min, eps_max = 1e-4, 1e-1
    cmap = 'inferno'
    
    def log_fmt(x, pos):
        if x <= 0: return ""
        return f"{np.log10(x):.0f}"
    formatter = FuncFormatter(log_fmt)
    
    axes_flat = axes.flatten()
    
    im = None
    for i, ax in enumerate(axes_flat):
        if i >= len(data_list):
            ax.set_visible(False)
            continue
            
        ax.tick_params(which='major', length=10, width=2.0, labelsize=18, color='black', pad=8)
        ax.tick_params(which='minor', length=5, width=1.2, color='black')
        
        row = i // cols
        col = i % cols
        
        if col == 0:
            ax.set_ylabel(r"$\log_{10}(r \,[\rm AU])$", fontsize=20, labelpad=10)
        
        if row == rows - 1:
            ax.set_xlabel(r"$\log_{10}(t \,[\rm Myr])$", fontsize=20, labelpad=10)
            
        data = data_list[i]
        title = titles[i]
        
        
        if row == 0 and col_titles is not None and col < len(col_titles):
            ax.set_title(col_titles[col], fontsize=22, pad=15)
            
        if title:
            ax.text(0.05, 0.95, title, transform=ax.transAxes, color='white', 
                    fontsize=14, fontweight='bold', va='top', ha='left',
                    bbox=dict(facecolor='black', alpha=0.6, edgecolor='none', pad=6))

        if data is None:
            ax.set_facecolor('black')
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.xaxis.set_major_formatter(formatter)
            ax.yaxis.set_major_formatter(formatter)
            ax.set_xlim([1e-3, 10.0])
            ax.set_ylim([0.7, 100.0])
            continue
            
        t_arr = data['t_arr'] / (3.15576e7 * 1e6)
        r_arr = data['r_arr']
        eps_mat = data['eps_mat']
        
        eps_mat_plot = np.where(eps_mat < eps_min, eps_min, eps_mat)
        im = ax.pcolormesh(t_arr, r_arr, eps_mat_plot.T, cmap=cmap, norm=LogNorm(vmin=eps_min, vmax=eps_max), shading='nearest')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.xaxis.set_major_formatter(formatter)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_xlim([1e-3, 10.0])
        ax.set_ylim([0.7, 100.0])
        
        AU_cm = 1.495978707e13
        r_ice_real = np.array([r_snow_time_cgs(t_sec) for t_sec in data['t_arr']]) / AU_cm
        ax.plot(t_arr, r_ice_real, color='cyan', linestyle=':', alpha=1, lw=2.5)
        
        # Remove first x tick for all columns except the first one
        if col > 0 and row == rows - 1:
            ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            ax.set_xticklabels(['', r'$-2$', r'$-1$', r'$0$', r'$1$'])

    if im is not None:
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.89, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(im, cax=cbar_ax, extend='both')
        cbar.set_label(r"Relación Polvo-Gas $\epsilon$", fontsize=26, rotation=270, labelpad=40)
        cbar.ax.tick_params(labelsize=22)
    
    if suptitle:
        fig.suptitle(suptitle, fontsize=24, y=0.96 if rows > 1 else 1.0)
    
    base_name = os.path.splitext(out_name)[0]
    out_png = os.path.join(DIRECTORIO_SALIDA, f"{base_name}.png")
    out_pdf = os.path.join(DIRECTORIO_SALIDA, f"{base_name}.pdf")
    plt.savefig(out_png, dpi=350, bbox_inches='tight')
    plt.savefig(out_pdf, dpi=350, bbox_inches='tight')
    plt.close()
    print(f" -> Guardado: {out_png} y {out_pdf}\n")

def figura1():
    print("=== Generando Figura 1: Comparación v_frag 1, 3, 10 m/s (3x3) ===")
    paths = [
        find_run_path("smooth", 1, "0.001"),
        find_run_path("smooth", 3, "0.001"),
        find_run_path("smooth", 10, "0.001"),
        find_run_path("gap", 1, "0.001", "1.0"),
        find_run_path("gap", 3, "0.001", "1.0"),
        find_run_path("gap", 10, "0.001", "1.0"),
        find_run_path("sinusoidal", 1, "0.001", n_gap=5, amp=1.0),
        find_run_path("sinusoidal", 3, "0.001", n_gap=5, amp=1.0),
        find_run_path("sinusoidal", 10, "0.001", n_gap=5, amp=1.0)
    ]
    col_titles = [
        r"$v_{\rm frag}=1$ m/s",
        r"$v_{\rm frag}=3$ m/s",
        r"$v_{\rm frag}=10$ m/s"
    ]
    titles = [None] * 9
    data_list = [extract_data(p) for p in paths]
    plot_hovmoller_grid(data_list, titles, "figura1_vfrags.png", rows=3, cols=3, figsize=(18, 12), suptitle=r"$\alpha=10^{-3}$", col_titles=col_titles)

def figura3():
    print("=== Generando Figura 3: Mosaico de Turbulencia (Uni Gap vs Sinusoidal) ===")
    alphas = ["0.0001", "0.0005", "0.001"]
    col_titles = [r"$\alpha=10^{-4}$", r"$\alpha=5\times10^{-4}$", r"$\alpha=10^{-3}$"]
    
    paths = []
    titles = []
    
    # Fila 1: Gap único 0.1 M_J
    for a in alphas:
        paths.append(find_run_path("gap", 10, a, "0.1"))
        titles.append(None)
        
    # Fila 2: Sinusoidal
    for a in alphas:
        paths.append(find_run_path("sinusoidal", 10, a, n_gap=5, amp=1.0))
        titles.append(None)
        
    data_list = [extract_data(p) for p in paths]
    plot_hovmoller_grid(data_list, titles, "figura3_turbulencia.png", rows=2, cols=3, figsize=(15, 8), suptitle=None, col_titles=col_titles)

def figura2():
    print("=== Generando Figura 2: Mosaico de Masas del Embrión ===")
    masas = ["0.01", "0.1", "0.5"]
    
    paths = []
    titles = []
    
    # Fila 1: Gap único
    for m in masas:
        paths.append(find_run_path("gap", 10, "0.001", m))
        titles.append(r"$M_{\rm gap}=" + str(m) + r" M_J$")
        
    # Fila 2: Sinusoidal (amplitudes)
    amps = ["0.5", "1.0", "3.0"]
    for a in amps:
        paths.append(find_run_path("sinusoidal", 10, "0.001", n_gap=5, amp=float(a)))
        titles.append(r"$A=" + str(a) + r"$")
        
    data_list = [extract_data(p) for p in paths]
    plot_hovmoller_grid(data_list, titles, "figura2_masas.png", rows=2, cols=3, figsize=(15, 8), suptitle=r"$\alpha=10^{-3}$")


def figura4_gap_unico():
    """
    Figura A1 — Gap unico 1x3 (v_frag=1,3,10 m/s).
    Colorbar inferno: M_gap [M_Jup]. Alpha=0.001, r_gap=10 AU.
    Salida: data/figures/paper/figura_A1_gap_unico.png
    """
    print("=== Generando Figura A1: Gap unico 1x3 ===")
    from plot_figura_vfrag_alpha import build_figura_A1
    build_figura_A1()


def figura5_sinusoidal():
    """
    Figura A2 — Sinusoidal 1x3 (v_frag=1,3,10 m/s).
    Colorbar viridis: Amplitud A. N_gaps=10, Alpha=0.001.
    Salida: data/figures/paper/figura_A2_sinusoidal.png
    """
    print("=== Generando Figura A2: Sinusoidal 1x3 ===")
    from plot_figura_vfrag_alpha import build_figura_A2
    build_figura_A2()


def figura6_regimenes_mgap():
    """
    Figura B — Regimenes de M_gap 2x2.
    Colorbar inferno: alpha. v_frag=10, r_gap=10 AU.
    Paneles: M_gap = 0.01, 0.1, 0.3, 0.5 M_Jup.
    Salida: data/figures/paper/figura_B_regimenes_mgap.png
    """
    print("=== Generando Figura B: Regimenes M_gap 2x2 ===")
    from plot_figura_vfrag_alpha import build_figura_B
    build_figura_B()


def figura7_A3_combinada():
    """
    Figura A3 — Gap unico + Sinusoidal combinadas 2x3 (v_frag=1,3,10 m/s).
    Fila 1: Gap unico (inferno: M_gap). Fila 2: Sinusoidal (viridis: A).
    Salida: data/figures/paper/figura_A3_gap_sinusoidal.png
    """
    print("=== Generando Figura A3: Gap unico + Sinusoidal 2x3 ===")
    from plot_figura_vfrag_alpha import build_figura_A3
    build_figura_A3()


def figura8_facet_combinado():
    """
    Figura C — Facet scatter combinado 2x4 (General + Sinusoidal).
    Columnas: alpha = 0.0001, 0.0003, 0.0005, 0.0007.
    Rangos: Earth analogs (verde) y Waterworlds (azul).
    """
    print("=== Generando Figura C: Facet Scatter combinado 2x4 (Lineal) ===")
    from plot_figura_vfrag_alpha import build_figura_C_facet
    build_figura_C_facet(x_scale='linear', x_max=1.05, suffix='_lineal')
    
    print("=== Generando Figura C: Facet Scatter combinado 2x4 (Log 0-10%) ===")
    build_figura_C_facet(x_scale='log', x_max=0.1, suffix='_log_10pct')


if __name__ == '__main__':
    print("\n==================================================")
    print("      GENERADOR DE FIGURAS PARA RESULTADOS")
    print("==================================================")
    
    figura1()
    figura2()
    figura3()
    figura4_gap_unico()
    figura5_sinusoidal()
    figura6_regimenes_mgap()
    figura7_A3_combinada()
    figura8_facet_combinado()
    
    print("\n!Todas las figuras han sido generadas con exito!")

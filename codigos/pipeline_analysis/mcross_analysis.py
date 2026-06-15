import os
import glob
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
import matplotlib.patheffects as pe

def generate_mcross_figures(runs_dir, out_dir):
    # Set publication style
    plt.style.use('default')
    plt.rcParams.update({
        'font.family': 'serif',
        'font.size': 16,
        'axes.labelsize': 18,
        'axes.titlesize': 18,
        'xtick.labelsize': 14,
        'ytick.labelsize': 14,
        'legend.fontsize': 14,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'xtick.minor.visible': True,
        'ytick.minor.visible': True,
    })

    vfrags = [1, 3, 10]
    all_data = []

    print("Recolectando datos masivos de todos los escenarios...")
    for vf in vfrags:
        vf_dir = os.path.join(runs_dir, f"vf_{vf}ms")
        if not os.path.isdir(vf_dir): continue
        
        # Iterar sobre todos los escenarios
        for s in os.listdir(vf_dir):
            if not os.path.isdir(os.path.join(vf_dir, s)): continue
            if s not in ['general', 'sinusoidal', 'delayed', 'rounded']: continue
            
            cache_files = glob.glob(os.path.join(vf_dir, s, "**/cache_*.pkl"), recursive=True)
            for cf in cache_files:
                try:
                    with open(cf, 'rb') as f:
                        data_list = pickle.load(f)
                        for d in data_list:
                            # Extract M_final
                            m_final = None
                            if 'm_emb' in d: m_final = d['m_emb'][-1]
                            elif 'mass_e' in d: m_final = d['mass_e'][-1]
                            elif 'M_final' in d: m_final = d['M_final']
                            
                            # Extract f_water
                            f_water = None
                            if 'f_water' in d: f_water = d['f_water'][-1]
                            elif 'frac_h2o_final' in d: f_water = d['frac_h2o_final'] / 100.0
                            elif 'frac_h2o_percent' in d: f_water = d['frac_h2o_percent'] / 100.0
                            
                            if m_final is None or f_water is None: continue
                            
                            # Extract M_cross
                            t = np.array(d.get('times_yr', []))
                            m = np.array(d.get('mass_e', []))
                            t_c = d.get('t_cross_1au', 0)
                            if t_c is None: t_c = 0
                            
                            if len(t) > 0 and len(m) > 0:
                                m_cross = m[t <= t_c][-1] if any(t <= t_c) else m[0]
                            else:
                                continue
                                
                            all_data.append({
                                'vfrag': vf,
                                'scenario': s,
                                'M_cross': m_cross,
                                'M_final': m_final,
                                'f_water': f_water
                            })
                except Exception:
                    pass

    df = pd.DataFrame(all_data)
    if df.empty:
        print("No se encontraron datos.")
        return

    os.makedirs(out_dir, exist_ok=True)
    
    # -------------------------------------------------------------
    # Figura C: M_final vs M_cross (Color: f_water)
    # -------------------------------------------------------------
    print("Generando Figura C...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0)
    
    cmap_c = plt.get_cmap('inferno')
    vmin_c, vmax_c = 0.0, df['f_water'].max()
    if pd.isna(vmax_c) or vmax_c <= 0: vmax_c = 1.0
    
    sc_c = None
    for i, vf in enumerate(vfrags):
        ax = axes[i]
        sub = df[df['vfrag'] == vf]
        if not sub.empty:
            sc_c = ax.scatter(sub['M_cross'], sub['M_final'], c=sub['f_water'], 
                              cmap=cmap_c, vmin=vmin_c, vmax=vmax_c, s=50, alpha=0.8, edgecolors='none', zorder=3)
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        ax.tick_params(which='major', direction='in', top=True, right=True, bottom=True, left=True, length=8, width=1.0)
        ax.tick_params(which='minor', direction='in', top=True, right=True, bottom=True, left=True, length=4, width=0.8)
        
        ax.set_title(rf"$v_{{frag}} = {vf}$ m/s", pad=10)
        ax.set_xlabel(r"Masa al cruzar la snowline ($M_{\rm cross}$) [$M_\oplus$]")
        if i == 0:
            ax.set_ylabel(r"Masa Final ($M_{\rm final}$) [$M_\oplus$]")
            
        ax.plot([1e-4, 10], [1e-4, 10], 'k--', alpha=0.3, zorder=1) # Línea 1:1
            
    if sc_c is not None:
        cbar = fig.colorbar(sc_c, ax=axes.ravel().tolist(), pad=0.02, aspect=30)
        cbar.set_label(r"Fracción de Agua Final ($f_{\rm H_2O}$)", fontsize=18)
        cbar.ax.tick_params(which='major', direction='in', length=8)
        cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
        cbar.ax.tick_params(which='minor', direction='in', length=4)
        
    fig.savefig(os.path.join(out_dir, "FigC_Mfinal_vs_Mcross.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("-> Guardado: FigC_Mfinal_vs_Mcross.png")
    
    # -------------------------------------------------------------
    # Figura D: f_water vs M_cross (Color: M_final)
    # -------------------------------------------------------------
    print("Generando Figura D...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0)
    
    cmap_d = plt.get_cmap('viridis')
    vmin_m = max(df['M_final'].min(), 1e-3)
    vmax_m = df['M_final'].max()
    norm = mcolors.LogNorm(vmin=vmin_m, vmax=vmax_m)
    
    sc_d = None
    for i, vf in enumerate(vfrags):
        ax = axes[i]
        sub = df[df['vfrag'] == vf]
        if not sub.empty:
            sc_d = ax.scatter(sub['M_cross'], sub['f_water'], c=sub['M_final'], 
                              cmap=cmap_d, norm=norm, s=50, alpha=0.8, edgecolors='none', zorder=3)
        
        ax.set_xscale('log')
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.tick_params(which='major', direction='in', top=True, right=True, bottom=True, left=True, length=8, width=1.0)
        ax.tick_params(which='minor', direction='in', top=True, right=True, bottom=True, left=True, length=4, width=0.8)
        
        ax.set_title(rf"$v_{{frag}} = {vf}$ m/s", pad=10)
        ax.set_xlabel(r"Masa al cruzar la snowline ($M_{\rm cross}$) [$M_\oplus$]")
        if i == 0:
            ax.set_ylabel(r"Fracción de Agua Final ($f_{\rm H_2O}$)")
            
    if sc_d is not None:
        cbar = fig.colorbar(sc_d, ax=axes.ravel().tolist(), pad=0.02, aspect=30)
        cbar.set_label(r"Masa Final ($M_{\rm final}$) [$M_\oplus$]", fontsize=18)
        cbar.ax.tick_params(which='major', direction='in', length=8)
        cbar.ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        cbar.ax.tick_params(which='minor', direction='in', length=4)
        
    fig.savefig(os.path.join(out_dir, "FigD_Agua_vs_Mcross.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("-> Guardado: FigD_Agua_vs_Mcross.png")

    # -------------------------------------------------------------
    # Figura E: f_water vs M_cross/M_final (Color: M_final)
    # -------------------------------------------------------------
    print("Generando Figura E...")
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    plt.subplots_adjust(wspace=0.0)
    
    sc_e = None
    for i, vf in enumerate(vfrags):
        ax = axes[i]
        sub = df[df['vfrag'] == vf]
        if not sub.empty:
            ratio = sub['M_cross'] / sub['M_final']
            sc_e = ax.scatter(ratio, sub['f_water'], c=sub['M_final'], 
                              cmap=cmap_d, norm=norm, s=50, alpha=0.8, edgecolors='none', zorder=3)
        
        ax.set_xscale('log')
        ax.minorticks_on()
        ax.xaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
        ax.tick_params(which='major', direction='in', top=True, right=True, bottom=True, left=True, length=8, width=1.0)
        ax.tick_params(which='minor', direction='in', top=True, right=True, bottom=True, left=True, length=4, width=0.8)
        
        ax.set_title(rf"$v_{{frag}} = {vf}$ m/s", pad=10)
        ax.set_xlabel(r"Fracción de Masa al Cruce ($M_{\rm cross} / M_{\rm final}$)")
        if i == 0:
            ax.set_ylabel(r"Fracción de Agua Final ($f_{\rm H_2O}$)")
            
    if sc_e is not None:
        cbar = fig.colorbar(sc_e, ax=axes.ravel().tolist(), pad=0.02, aspect=30)
        cbar.set_label(r"Masa Final ($M_{\rm final}$) [$M_\oplus$]", fontsize=18)
        cbar.ax.tick_params(which='major', direction='in', length=8)
        cbar.ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10)*1.0, numticks=15))
        cbar.ax.tick_params(which='minor', direction='in', length=4)
        
    fig.savefig(os.path.join(out_dir, "FigE_Agua_vs_Ratio.png"), dpi=300, bbox_inches='tight')
    plt.close(fig)
    print("-> Guardado: FigE_Agua_vs_Ratio.png")

if __name__ == "__main__":
    runs_dir = os.path.join(os.path.dirname(__file__), "..", "data", "runs")
    out_dir = os.path.join(os.path.dirname(__file__), "..", "data", "figures", "benchmarks")
    generate_mcross_figures(runs_dir, out_dir)

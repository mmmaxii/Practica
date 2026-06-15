import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from matplotlib.lines import Line2D

# Configuración global de fuentes
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.titlesize': 18
})

class PlotterGeneral:
    @staticmethod
    def plot_lines_grouped_by_rgap(data, fig_dir, alpha_val=None):
        groups = {}
        for item in data:
            rg = item['r_gap']
            if rg not in groups: groups[rg] = []
            groups[rg].append(item)
            
        r_gaps = sorted(list(groups.keys()))
        if not r_gaps: return
        
        fig, axes = plt.subplots(2, 2, figsize=(14, 11))
        title = "Evolución Temporal"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        fig.suptitle(title, fontsize=18)
        
        axes = axes.flatten()
        cmap = plt.get_cmap('viridis')
        all_mgaps = sorted(list(set([item['M_gap'] for item in data])))
        norm = LogNorm(vmin=min(all_mgaps), vmax=max(all_mgaps))
        
        for i, (ax, rg) in enumerate(zip(axes, r_gaps)):
            ax.set_title(f"Gap a {rg} AU")
            runs_rg = sorted(groups[rg], key=lambda x: x['M_gap'])
            t_cross_avg = []
            
            for item in runs_rg:
                color = cmap(norm(item['M_gap']))
                ls = '-' if item.get('completed', True) else '-.'
                
                t_myr_log = np.log10(item['times_yr'] / 1e6)
                mass_e_log = np.log10(item['mass_e'])
                
                # M_iso true curve from PA3Py
                if 'm_iso_e' in item:
                    m_iso_log = np.log10(item['m_iso_e'])
                    ax.plot(t_myr_log, m_iso_log, color=color, ls=':', alpha=0.4, lw=1.5, zorder=1)
                
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
                custom_lines = [Line2D([0], [0], color=cmap(norm(mg)), lw=2.5) for mg in all_mgaps]
                labels = [f"{mg} $M_{{\mathrm{{Jup}}}}$" for mg in all_mgaps]
                if any(not item.get('completed', True) for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls='-.'))
                    labels.append("Incompleta (< 10Myr)")
                if any('m_iso_e' in item for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=1.5, ls=':', alpha=0.5))
                    labels.append("M_iso (límite)")
                if t_cross_avg:
                    custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
                    labels.append("Snowline a 1au")
                ax.legend(custom_lines, labels, loc="upper left", fontsize=11)

        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_rgap.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_lines_grouped_by_mgap(data, fig_dir, alpha_val=None):
        groups = {}
        for item in data:
            mg = item['M_gap']
            if mg not in groups: groups[mg] = []
            groups[mg].append(item)
            
        m_gaps = sorted(list(groups.keys()))
        if not m_gaps: return
        n_mg = len(m_gaps)
        cols = 3
        rows = (n_mg + cols - 1) // cols
        
        fig, axes = plt.subplots(rows, cols, figsize=(16, 5 * rows))
        title = "Evolución Temporal por Profundidad"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        fig.suptitle(title, fontsize=18)
        
        if not isinstance(axes, np.ndarray): axes = np.array([axes])
        axes = axes.flatten()
        cmap = plt.get_cmap('plasma')
        all_rgaps = sorted(list(set([item['r_gap'] for item in data])))
        norm = Normalize(vmin=min(all_rgaps), vmax=max(all_rgaps))
        
        for i, (ax, mg) in enumerate(zip(axes, m_gaps)):
            ax.set_title(f"Profundidad: {mg} $M_{{\mathrm{{Jup}}}}$")
            runs_mg = sorted(groups[mg], key=lambda x: x['r_gap'])
            t_cross_avg = []
            
            for k, item in enumerate(runs_mg):
                color = cmap(norm(item['r_gap']))
                completed = item.get('completed', True)
                ls = '-' if completed else '-.'
                
                t_myr_log = np.log10(item['times_yr'] / 1e6)
                mass_e_log = np.log10(item['mass_e'])
                
                if 'm_iso_e' in item:
                    m_iso_log = np.log10(item['m_iso_e'])
                    ax.plot(t_myr_log, m_iso_log, color=color, ls=':', alpha=0.4, lw=1.5, zorder=1)
                
                mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
                
                lw = 2.0 + k * 2.0 if mg == 0.01 else 2.5
                z_ord = 10 - k if mg == 0.01 else 2
                
                ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=lw, zorder=z_ord, ls=ls)
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
                ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0)
                
            row, col = i // cols, i % cols
            if col == 0: ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
            last_in_col = max([idx for idx in range(n_mg) if idx % cols == col])
            if i == last_in_col: ax.set_xlabel("Log(t [Myr])")
            
            if i == 0:
                custom_lines = [Line2D([0], [0], color=cmap(norm(rg)), lw=2.5) for rg in all_rgaps]
                labels = [f"{rg} AU" for rg in all_rgaps]
                if any(not item.get('completed', True) for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls='-.'))
                    labels.append("Incompleta (< 10Myr)")
                if any('m_iso_e' in item for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=1.5, ls=':', alpha=0.5))
                    labels.append("M_iso (límite)")
                ax.legend(custom_lines, labels, loc="upper left", fontsize=11)
                
        for j in range(n_mg, len(axes)): axes[j].set_visible(False)
        plt.tight_layout()
        plt.subplots_adjust(top=0.92)
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_mgap.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_heatmaps(data, fig_dir, alpha_val=None):
        r_gaps = sorted(list(set([d['r_gap'] for d in data])))
        m_gaps = sorted(list(set([d['M_gap'] for d in data])))
        
        Z_mass = np.full((len(m_gaps), len(r_gaps)), np.nan)
        Z_h2o = np.full((len(m_gaps), len(r_gaps)), np.nan)
        
        for item in data:
            j, i = r_gaps.index(item['r_gap']), m_gaps.index(item['M_gap'])
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(r_gaps)), np.arange(len(m_gaps))
        
        # Escalar el tamaño de la figura para que no se apriete el texto
        width = max(9, len(r_gaps) * 1.3)
        height = max(7, len(m_gaps) * 1.1)
        
        # Heatmap Masa
        fig, ax = plt.subplots(figsize=(width, height))
        vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
        ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
        for i in range(len(m_gaps)):
            for j in range(len(r_gaps)):
                val = Z_mass[i, j]
                if not np.isnan(val):
                    color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                    ax.text(j, i, rf"{val:.3f} $M_\oplus$", ha="center", va="center", color=color, fontsize=10)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(r_gaps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(m_gaps)
        ax.set_xlabel("Posición del Gap [AU]")
        ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
        title_mass = "Masa Final del Embrión"
        if alpha_val is not None: title_mass += rf" ($\alpha = {alpha_val}$)"
        ax.set_title(title_mass)
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_masa.png"), dpi=200, bbox_inches='tight')
        plt.close()
        
        # Heatmap H2O
        fig, ax = plt.subplots(figsize=(width, height))
        vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
        if vmax_h2o == 0: vmax_h2o = 1.0
        ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', vmin=0, vmax=vmax_h2o)
        for i in range(len(m_gaps)):
            for j in range(len(r_gaps)):
                val = Z_h2o[i, j]
                if not np.isnan(val):
                    color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                    ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=color, fontsize=11)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(r_gaps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(m_gaps)
        ax.set_xlabel("Posición del Gap [AU]")
        ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]")
        title_h2o = "Fracción de Agua Final"
        if alpha_val is not None: title_h2o += rf" ($\alpha = {alpha_val}$)"
        ax.set_title(title_h2o)
        plt.savefig(os.path.join(fig_dir, "heatmap_h2o.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_categorical_map(data, fig_dir, alpha_val=None):
        import matplotlib.patches as mpatches
        
        r_gaps = sorted(list(set([d['r_gap'] for d in data])))
        m_gaps = sorted(list(set([d['M_gap'] for d in data])))
        
        Z_mass = np.full((len(m_gaps), len(r_gaps)), np.nan)
        Z_h2o = np.full((len(m_gaps), len(r_gaps)), np.nan)
        
        for item in data:
            j, i = r_gaps.index(item['r_gap']), m_gaps.index(item['M_gap'])
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(r_gaps)), np.arange(len(m_gaps))
        
        # Escalar el tamaño de la figura para que no se apriete el texto
        width = max(10, len(r_gaps) * 1.3)
        height = max(8, len(m_gaps) * 1.1)
        fig, ax = plt.subplots(figsize=(width, height))
        
        # Fondo gris claro para "Fracaso"
        ax.set_facecolor('#e0e0e0')
        
        # Enmascarar donde la masa es menor a 0.1
        Z_h2o_masked = np.ma.masked_where(Z_mass < 0.1, Z_h2o)
        
        # Plot continuous colormap for successful embryos
        cmap = plt.get_cmap('PuBu')
        # Anclar vmax a 50.0 universalmente para comparar entre mapas de distintos alphas
        vmax_h2o_universal = 50.0
        
        im = ax.imshow(Z_h2o_masked, origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=vmax_h2o_universal)
        
        # Add colorbar
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label(r"Fracción de Agua Final (%)", fontsize=14)
        
        # Draw exact pixel boundaries separating success from failure (step-like contour)
        for i in range(len(m_gaps)):
            for j in range(len(r_gaps)):
                mass = Z_mass[i, j]
                if not np.isnan(mass) and mass >= 0.1:
                    # Check bottom edge
                    if i == 0 or np.isnan(Z_mass[i-1, j]) or Z_mass[i-1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i-0.5, i-0.5], color='black', lw=3)
                    # Check top edge
                    if i == len(m_gaps)-1 or np.isnan(Z_mass[i+1, j]) or Z_mass[i+1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i+0.5, i+0.5], color='black', lw=3)
                    # Check left edge
                    if j == 0 or np.isnan(Z_mass[i, j-1]) or Z_mass[i, j-1] < 0.1:
                        ax.plot([j-0.5, j-0.5], [i-0.5, i+0.5], color='black', lw=3)
                    # Check right edge
                    if j == len(r_gaps)-1 or np.isnan(Z_mass[i, j+1]) or Z_mass[i, j+1] < 0.1:
                        ax.plot([j+0.5, j+0.5], [i-0.5, i+0.5], color='black', lw=3)
        
        # Add text to each cell
        for i in range(len(m_gaps)):
            for j in range(len(r_gaps)):
                mass = Z_mass[i, j]
                water = Z_h2o[i, j]
                if not np.isnan(mass) and not np.isnan(water):
                    # Solo imprimimos texto en el "Continente del Éxito" (Limpia el cementerio)
                    if mass >= 0.1:
                        text_color = 'white' if water > (vmax_h2o_universal * 0.5) else 'black'
                        # Elimina la redundancia (sin M_Tierra ni H2O)
                        label = f"{mass:.2f}\n{water:.1f}%"
                        ax.text(j, i, label, ha="center", va="center", color=text_color, fontsize=11, fontweight='bold')
                    
        ax.set_xticks(x_idx)
        ax.set_xticklabels(r_gaps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(m_gaps)
        ax.set_xlabel("Posición del Gap [AU]", fontsize=14)
        ax.set_ylabel(r"Profundidad del Gap [$M_{\mathrm{Jup}}$]", fontsize=14)
        ax.tick_params(labelsize=12)
        
        # Add explanation of units at the bottom
        ax.annotate("Valores en celdas: Masa final [$M_\\oplus$] (arriba) | Fracción de Agua (abajo)", 
                    xy=(0.5, -0.15), xycoords='axes fraction', ha='center', va='top', fontsize=12)
        
        title = "Mapa de Evolución (Continente del Éxito)"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        ax.set_title(title, fontsize=16)
        
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_categorico.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_lines_grouped_by_rgap_mosaic(data, fig_dir, alpha_val=None):
        import math
        groups = {}
        for item in data:
            rg = item['r_gap']
            if rg not in groups: groups[rg] = []
            groups[rg].append(item)
            
        r_gaps = sorted(list(groups.keys()))
        if not r_gaps: return
        n_rg = len(r_gaps)
        cols = min(n_rg, 3)
        rows = math.ceil(n_rg / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.90)
        
        if not isinstance(axes, np.ndarray): axes = np.array([axes])
        axes = axes.flatten()
        cmap = plt.get_cmap('viridis')
        all_mgaps = sorted(list(set([item['M_gap'] for item in data])))
        norm = LogNorm(vmin=min(all_mgaps), vmax=max(all_mgaps))
        
        min_m = np.nanmin([np.nanmin(item['mass_e']) for item in data])
        if min_m > 1e-2: min_m = 1e-2
        min_m_plot = min_m / 3.0
        max_m_plot = 10.0
        
        for i, (ax, rg) in enumerate(zip(axes, r_gaps)):
            ax.text(0.05, 0.95, r"$r_{\rm gap} = " + f"{rg}" + r"\,{\rm AU}$", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=14,
                    bbox=dict(facecolor='white', edgecolor='lightgray', alpha=0.8))
            
            ax.tick_params(direction='in', top=True, right=True, which='major', labelsize=14, length=6, width=1.5)
            ax.tick_params(direction='in', top=True, right=True, which='minor', length=3, width=0.5)
            runs_rg = sorted(groups[rg], key=lambda x: x['M_gap'])
            t_cross_avg = []
            
            styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
            for item in runs_rg:
                color = cmap(norm(item['M_gap']))
                idx = all_mgaps.index(item['M_gap'])
                ls = styles_list[idx % len(styles_list)]
                t_myr = item['times_yr'] / 1e6
                mass_e = item['mass_e']
                if 'm_iso_e' in item:
                    ax.plot(t_myr, item['m_iso_e'], color='saddlebrown', ls='--', alpha=0.7, lw=2.5, zorder=1)
                
                mask = np.isfinite(t_myr) & np.isfinite(mass_e)
                t_valid = t_myr[mask]
                m_valid = mass_e[mask]
                if len(t_valid) > 0 and t_valid[-1] < 10.0:
                    t_valid = np.append(t_valid, 10.0)
                    m_valid = np.append(m_valid, m_valid[-1])
                    
                ax.plot(t_valid, m_valid, color=color, alpha=0.8, lw=3.5, ls=ls, zorder=2)
                if item['t_cross_1au']: t_cross_avg.append(item['t_cross_1au'] / 1e6)
                
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([1e-3, 10.0])
            
            import matplotlib.ticker as ticker
            row, col = i // cols, i % cols
            if col == 0:
                ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            else:
                ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
            ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
            
            ax.set_ylim([min_m_plot, max_m_plot])
            ax.grid(True, which='major', color='lightgray', linestyle='--', alpha=0.5)
            ax.grid(True, which='minor', color='whitesmoke', linestyle=':', alpha=0.3)
            ax.axhline(0.1, color='black', ls=':', alpha=0.8, lw=2.0)
            
            if t_cross_avg:
                t_med = np.median(t_cross_avg)
                ax.axvline(t_med, color='blue', ls='--', alpha=0.7, lw=2.0)
                
            row, col = i // cols, i % cols
            if col == 0: ax.set_ylabel(r"Masa del embrión [$M_\oplus$]", fontsize=16)
            last_in_col = max([idx for idx in range(n_rg) if idx % cols == col])
            if i == last_in_col or i + cols >= len(axes): ax.set_xlabel("Tiempo [Myr]", fontsize=16)
            
        for j in range(n_rg, len(axes)):
            axes[j].set_visible(False)
            
        styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
        custom_lines = [Line2D([0], [0], color=cmap(norm(mg)), lw=3.5, ls=styles_list[i % len(styles_list)]) for i, mg in enumerate(all_mgaps)]
        labels = [f"{mg} $M_{{\mathrm{{Jup}}}}$" for mg in all_mgaps]
        custom_lines.append(Line2D([0], [0], color='black', lw=2.0, ls=':', alpha=0.8))
        labels.append("Umbral (0.1 $M_\oplus$)")
        if any('m_iso_e' in item for item in data):
            custom_lines.append(Line2D([0], [0], color='saddlebrown', lw=2.5, ls='--', alpha=0.7))
            labels.append("$M_{\mathrm{iso}}$")
            
        ncols = min(len(labels), 6)
        fig.legend(custom_lines, labels, loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=ncols, frameon=False, fontsize=14)
            
        title = "Evolución Temporal por Posición (Mosaico)"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        fig.suptitle(title, fontsize=16, y=0.98)
        
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_rgap_mosaic.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_lines_grouped_by_mgap_mosaic(data, fig_dir, alpha_val=None):
        import math
        groups = {}
        for item in data:
            mg = item['M_gap']
            if mg not in groups: groups[mg] = []
            groups[mg].append(item)
            
        m_gaps = sorted(list(groups.keys()))
        if not m_gaps: return
        n_mg = len(m_gaps)
        cols = min(n_mg, 3)
        rows = math.ceil(n_mg / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        
        if not isinstance(axes, np.ndarray): axes = np.array([axes])
        axes = axes.flatten()
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list("custom_dist", ["darkred", "darkorange", "steelblue", "navy"])
        all_rgaps = sorted(list(set([item['r_gap'] for item in data])))
        from matplotlib.colors import Normalize
        norm = Normalize(vmin=min(all_rgaps), vmax=max(all_rgaps))
        
        min_m = np.nanmin([np.nanmin(item['mass_e']) for item in data])
        if min_m > 1e-2: min_m = 1e-2
        min_m_plot = min_m / 3.0
        max_m_plot = 10.0
        
        for i, (ax, mg) in enumerate(zip(axes, m_gaps)):
            ax.text(0.05, 0.95, r"$M_{\rm gap} = " + f"{mg}" + r"\,M_{\rm Jup}$", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=14,
                    bbox=dict(facecolor='white', edgecolor='lightgray', alpha=0.8))
            
            ax.tick_params(direction='in', top=True, right=True, which='major', labelsize=14, length=6, width=1.5)
            ax.tick_params(direction='in', top=True, right=True, which='minor', length=3, width=0.5)
            runs_mg = sorted(groups[mg], key=lambda x: x['r_gap'])
            t_cross_avg = []
            
            styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
            for item in runs_mg:
                color = cmap(norm(item['r_gap']))
                idx = all_rgaps.index(item['r_gap'])
                ls = styles_list[idx % len(styles_list)]
                t_myr = item['times_yr'] / 1e6
                mass_e = item['mass_e']
                if 'm_iso_e' in item:
                    ax.plot(t_myr, item['m_iso_e'], color='saddlebrown', ls='--', alpha=0.7, lw=2.0, zorder=1)
                
                mask = np.isfinite(t_myr) & np.isfinite(mass_e)
                t_valid = t_myr[mask]
                m_valid = mass_e[mask]
                if len(t_valid) > 0 and t_valid[-1] < 10.0:
                    t_valid = np.append(t_valid, 10.0)
                    m_valid = np.append(m_valid, m_valid[-1])
                    
                ax.plot(t_valid, m_valid, color=color, alpha=0.8, lw=3.5, zorder=2, ls=ls)
                if item['t_cross_1au']: t_cross_avg.append(item['t_cross_1au'] / 1e6)
                
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([1e-3, 10.0])
            
            import matplotlib.ticker as ticker
            row, col = i // cols, i % cols
            if col == 0:
                ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            else:
                ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
            ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
            
            ax.set_ylim([min_m_plot, max_m_plot])
            ax.grid(True, which='major', color='lightgray', linestyle='--', alpha=0.5)
            ax.grid(True, which='minor', color='whitesmoke', linestyle=':', alpha=0.3)
            ax.axhline(0.1, color='black', ls=':', alpha=0.8, lw=2.0)
            
            if t_cross_avg:
                t_med = np.median(t_cross_avg)
                ax.axvline(t_med, color='blue', ls='--', alpha=0.7, lw=2.0)
                
            row, col = i // cols, i % cols
            if col == 0: ax.set_ylabel(r"Masa del embrión [$M_\oplus$]", fontsize=16)
            last_in_col = max([idx for idx in range(n_mg) if idx % cols == col])
            if i == last_in_col or i + cols >= len(axes): ax.set_xlabel("Tiempo [Myr]", fontsize=16)
            
        for j in range(n_mg, len(axes)):
            axes[j].set_visible(False)
            
        styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
        custom_lines = [Line2D([0], [0], color=cmap(norm(rg)), lw=3.5, ls=styles_list[i % len(styles_list)]) for i, rg in enumerate(all_rgaps)]
        labels = [f"{rg} AU" for rg in all_rgaps]
        custom_lines.append(Line2D([0], [0], color='black', lw=2.0, ls=':', alpha=0.8))
        labels.append("Umbral (0.1 $M_\oplus$)")
        if any('m_iso_e' in item for item in data):
            custom_lines.append(Line2D([0], [0], color='saddlebrown', lw=2.5, ls='--', alpha=0.7))
            labels.append("$M_{\mathrm{iso}}$")
            
        ncols = min(len(labels), 6)
        fig.legend(custom_lines, labels, loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=ncols, frameon=False, fontsize=14)
            
        title = "Evolución Temporal por Profundidad (Mosaico)"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        fig.suptitle(title, fontsize=16, y=0.98)
        
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_mgap_mosaic.png"), dpi=200, bbox_inches='tight')
        plt.close()

class PlotterSinusoidal:
    GAP_NAMES = {
        "PocosGaps": "Pocos Gaps (tipo HL Tau)",
        "Intermedio": "Intermedio",
        "MuchosGaps": "Muchos Gaps (tipo AS 209)"
    }
    GAP_ORDER = ["PocosGaps", "Intermedio", "MuchosGaps"]

    @staticmethod
    def plot_lines_grouped_by_gap_type(data, fig_dir):
        if all('n_gaps' in d for d in data):
            gtypes = sorted(list(set([d['n_gaps'] for d in data])))
            group_key = 'n_gaps'
            gap_names = {g: f"{g} Gaps" for g in gtypes}
            if 3 in gtypes: gap_names[3] = "3 Gaps (HL Tau)"
            if 5 in gtypes: gap_names[5] = "5 Gaps (Intermedio)"
            if 10 in gtypes: gap_names[10] = "10 Gaps"
            if 15 in gtypes: gap_names[15] = "15 Gaps"
            if 20 in gtypes: gap_names[20] = "20 Gaps (AS 209)"
        else:
            gtypes = PlotterSinusoidal.GAP_ORDER
            group_key = 'gap_type'
            gap_names = PlotterSinusoidal.GAP_NAMES
            
        groups = {g: [] for g in gtypes}
        for item in data:
            val = item.get(group_key)
            if val in groups:
                groups[val].append(item)
                
        all_amps = sorted(list(set([item['amp_val'] for item in data])))
        if not all_amps: return
        
        n_cols = len(gtypes)
        fig, axes = plt.subplots(1, n_cols, figsize=(6 * n_cols, 6), sharey=True)
        if n_cols == 1: axes = [axes]
        fig.suptitle("Evolución de Masa - Discos Sinusoidales", fontsize=20, y=1.05)
        
        cmap = plt.get_cmap('plasma')
        norm = LogNorm(vmin=min(all_amps), vmax=max(all_amps))
        
        for i, gtype in enumerate(gtypes):
            ax = axes[i]
            ax.set_title(gap_names.get(gtype, gtype), pad=10)
            runs_gt = sorted(groups[gtype], key=lambda x: x['amp_val'])
            t_cross_avg = []
            
            for item in runs_gt:
                color = cmap(norm(item['amp_val']))
                completed = item.get('completed', True)
                ls = '-' if completed else '-.'
                
                t_myr_log = np.log10(item['times_yr'] / 1e6)
                mass_e_log = np.log10(item['mass_e'])
                
                if 'm_iso_e' in item:
                    m_iso_log = np.log10(item['m_iso_e'])
                    ax.plot(t_myr_log, m_iso_log, color=color, ls=':', alpha=0.4, lw=1.5, zorder=1)
                
                mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
                ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=2.5, ls=ls)
                if item['t_cross_1au']: t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
                    
            ax.set_xlim([-3, 1])
            ax.grid(True, alpha=0.3)
            ax.set_xlabel("Log(t [Myr])")
            
            if t_cross_avg:
                t_med = np.median(t_cross_avg)
                ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0)
                
            if i == 0:
                ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
                custom_lines = [Line2D([0], [0], color=cmap(norm(amp)), lw=2.5) for amp in all_amps]
                labels = [f"A = {amp}" for amp in all_amps]
                if any(not item.get('completed', True) for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=2.5, ls='-.'))
                    labels.append("Incompleta (< 10Myr)")
                if any('m_iso_e' in item for item in data):
                    custom_lines.append(Line2D([0], [0], color='gray', lw=1.5, ls=':', alpha=0.5))
                    labels.append("M_iso (límite)")
                ax.legend(custom_lines, labels, loc="upper left", fontsize=10)
                
        min_m_global = np.nanmin([np.nanmin(np.log10(item['mass_e'])) for item in data])
        max_m_global = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in data])
        padding = (max_m_global - min_m_global) * 0.1
        axes[0].set_ylim([min_m_global - padding, max_m_global + padding])

        plt.tight_layout()
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_sinusoidal.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_heatmaps(data, fig_dir):
        if all('n_gaps' in d for d in data):
            gtypes = sorted(list(set([d['n_gaps'] for d in data])))
            group_key = 'n_gaps'
            gap_names = {g: f"{g} Gaps" for g in gtypes}
            if 3 in gtypes: gap_names[3] = "3 Gaps (HL Tau)"
            if 5 in gtypes: gap_names[5] = "5 Gaps (Intermedio)"
            if 10 in gtypes: gap_names[10] = "10 Gaps"
            if 15 in gtypes: gap_names[15] = "15 Gaps"
            if 20 in gtypes: gap_names[20] = "20 Gaps (AS 209)"
        else:
            gtypes = PlotterSinusoidal.GAP_ORDER
            group_key = 'gap_type'
            gap_names = PlotterSinusoidal.GAP_NAMES

        all_amps = sorted(list(set([d['amp_val'] for d in data])))
        Z_mass = np.full((len(gtypes), len(all_amps)), np.nan)
        Z_h2o = np.full((len(gtypes), len(all_amps)), np.nan)
        
        for item in data:
            val = item.get(group_key)
            if val not in gtypes: continue
            i, j = gtypes.index(val), all_amps.index(item['amp_val'])
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(all_amps)), np.arange(len(gtypes))
        y_labels = [gap_names.get(g, g) for g in gtypes]
        
        # Heatmap Masa
        height = max(4, len(gtypes) * 0.8)
        fig, ax = plt.subplots(figsize=(10, height))
        vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
        ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
        for i in range(len(gtypes)):
            for j in range(len(all_amps)):
                val = Z_mass[i, j]
                if not np.isnan(val):
                    color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                    ax.text(j, i, rf"{val:.2f} $M_\oplus$", ha="center", va="center", color=color, fontsize=11)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(all_amps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(y_labels)
        ax.set_xlabel("Amplitud de la perturbación (A)")
        ax.set_title("Masa Final del Embrión - Discos Sinusoidales")
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_masa.png"), dpi=200, bbox_inches='tight')
        plt.close()
        
        # Heatmap H2O
        fig, ax = plt.subplots(figsize=(10, height))
        vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
        if vmax_h2o == 0: vmax_h2o = 1.0 
        ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', vmin=0, vmax=vmax_h2o)
        for i in range(len(gtypes)):
            for j in range(len(all_amps)):
                val = Z_h2o[i, j]
                if not np.isnan(val):
                    color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                    ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=color, fontsize=11)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(all_amps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(y_labels)
        ax.set_xlabel("Amplitud de la perturbación (A)")
        plt.savefig(os.path.join(fig_dir, "heatmap_h2o.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_categorical_map(data, fig_dir):
        if all('n_gaps' in d for d in data):
            gtypes = sorted(list(set([d['n_gaps'] for d in data])))
            group_key = 'n_gaps'
            gap_names = {g: f"{g} Gaps" for g in gtypes}
            if 3 in gtypes: gap_names[3] = "3 Gaps (HL Tau)"
            if 5 in gtypes: gap_names[5] = "5 Gaps (Intermedio)"
            if 10 in gtypes: gap_names[10] = "10 Gaps"
            if 15 in gtypes: gap_names[15] = "15 Gaps"
            if 20 in gtypes: gap_names[20] = "20 Gaps (AS 209)"
        else:
            gtypes = PlotterSinusoidal.GAP_ORDER
            group_key = 'gap_type'
            gap_names = PlotterSinusoidal.GAP_NAMES

        all_amps = sorted(list(set([d['amp_val'] for d in data])))
        Z_mass = np.full((len(gtypes), len(all_amps)), np.nan)
        Z_h2o = np.full((len(gtypes), len(all_amps)), np.nan)
        
        for item in data:
            val = item.get(group_key)
            if val not in gtypes: continue
            i, j = gtypes.index(val), all_amps.index(item['amp_val'])
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(all_amps)), np.arange(len(gtypes))
        y_labels = [gap_names.get(g, g) for g in gtypes]
        
        height = max(5, len(gtypes) * 0.9)
        fig, ax = plt.subplots(figsize=(10, height))
        
        ax.set_facecolor('#e0e0e0')
        Z_h2o_masked = np.ma.masked_where(Z_mass < 0.1, Z_h2o)
        
        cmap = plt.get_cmap('PuBu')
        vmax_h2o_universal = 50.0
        
        im = ax.imshow(Z_h2o_masked, origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=vmax_h2o_universal)
        
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label(r"Fracción de Agua Final (%)", fontsize=14)
        
        for i in range(len(gtypes)):
            for j in range(len(all_amps)):
                mass = Z_mass[i, j]
                if not np.isnan(mass) and mass >= 0.1:
                    if i == 0 or np.isnan(Z_mass[i-1, j]) or Z_mass[i-1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i-0.5, i-0.5], color='black', lw=3)
                    if i == len(gtypes)-1 or np.isnan(Z_mass[i+1, j]) or Z_mass[i+1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i+0.5, i+0.5], color='black', lw=3)
                    if j == 0 or np.isnan(Z_mass[i, j-1]) or Z_mass[i, j-1] < 0.1:
                        ax.plot([j-0.5, j-0.5], [i-0.5, i+0.5], color='black', lw=3)
                    if j == len(all_amps)-1 or np.isnan(Z_mass[i, j+1]) or Z_mass[i, j+1] < 0.1:
                        ax.plot([j+0.5, j+0.5], [i-0.5, i+0.5], color='black', lw=3)
                        
        for i in range(len(gtypes)):
            for j in range(len(all_amps)):
                mass = Z_mass[i, j]
                water = Z_h2o[i, j]
                if not np.isnan(mass) and not np.isnan(water):
                    if mass >= 0.1:
                        text_color = 'white' if water > (vmax_h2o_universal * 0.5) else 'black'
                        label = f"{mass:.2f}\n{water:.1f}%"
                        ax.text(j, i, label, ha="center", va="center", color=text_color, fontsize=11, fontweight='bold')
                        
        ax.set_xticks(x_idx)
        ax.set_xticklabels(all_amps)
        ax.set_yticks(y_idx)
        ax.set_yticklabels(y_labels)
        ax.set_xlabel("Amplitud de la perturbación (A)", fontsize=14)
        ax.set_ylabel(r"$N_{\mathrm{gaps}}$", fontsize=14)
        ax.tick_params(labelsize=12)
        
        ax.annotate("Valores en celdas: Masa final [$M_\\oplus$] (arriba) | Fracción de Agua (abajo)", 
                    xy=(0.5, -0.15), xycoords='axes fraction', ha='center', va='top', fontsize=12)
        
        ax.set_title("Mapa de Evolución Sinusoidal (Continente del Éxito)", fontsize=16)
        
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_categorico.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_lines_grouped_by_amp_mosaic(data, fig_dir):
        import math
        groups = {}
        for item in data:
            amp = item['amp_val']
            if amp not in groups: groups[amp] = []
            groups[amp].append(item)
            
        amps = sorted(list(groups.keys()))
        if not amps: return
        n_amp = len(amps)
        cols = min(n_amp, 3)
        rows = math.ceil(n_amp / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.90)
        
        if not isinstance(axes, np.ndarray): axes = np.array([axes])
        axes = axes.flatten()
        cmap = plt.get_cmap('viridis')
        
        if all('n_gaps' in d for d in data):
            all_ngaps = sorted(list(set([item['n_gaps'] for item in data])))
            group_key = 'n_gaps'
        else:
            all_ngaps = PlotterSinusoidal.GAP_ORDER
            group_key = 'gap_type'
            
        from matplotlib.colors import LogNorm, Normalize
        if isinstance(all_ngaps[0], (int, float)):
            norm = LogNorm(vmin=min(all_ngaps), vmax=max(all_ngaps))
        else:
            norm = Normalize(vmin=0, vmax=max(1, len(all_ngaps)-1))
            
        min_m = np.nanmin([np.nanmin(item['mass_e']) for item in data])
        if min_m > 1e-2: min_m = 1e-2
        min_m_plot = min_m / 3.0
        max_m_plot = 10.0
        
        for i, (ax, amp) in enumerate(zip(axes, amps)):
            ax.text(0.05, 0.95, r"Amplitud = " + f"{amp}", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=14,
                    bbox=dict(facecolor='white', edgecolor='lightgray', alpha=0.8))
            
            ax.tick_params(direction='in', top=True, right=True, which='major', labelsize=14, length=6, width=1.5)
            ax.tick_params(direction='in', top=True, right=True, which='minor', length=3, width=0.5)
            
            runs_amp = sorted(groups[amp], key=lambda x: x[group_key])
            t_cross_avg = []
            
            styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
            for item in runs_amp:
                val = item[group_key]
                idx = all_ngaps.index(val)
                color = cmap(norm(val) if isinstance(val, (int, float)) else norm(idx))
                ls = styles_list[idx % len(styles_list)]
                
                t_myr = item['times_yr'] / 1e6
                mass_e = item['mass_e']
                if 'm_iso_e' in item:
                    ax.plot(t_myr, item['m_iso_e'], color='saddlebrown', ls='--', alpha=0.7, lw=2.5, zorder=1)
                
                mask = np.isfinite(t_myr) & np.isfinite(mass_e)
                t_valid = t_myr[mask]
                m_valid = mass_e[mask]
                if len(t_valid) > 0 and t_valid[-1] < 10.0:
                    t_valid = np.append(t_valid, 10.0)
                    m_valid = np.append(m_valid, m_valid[-1])
                    
                ax.plot(t_valid, m_valid, color=color, alpha=0.8, lw=3.5, ls=ls, zorder=2)
                if item['t_cross_1au']: t_cross_avg.append(item['t_cross_1au'] / 1e6)
                
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([1e-3, 10.0])
            
            import matplotlib.ticker as ticker
            row, col = i // cols, i % cols
            if col == 0:
                ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            else:
                ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
            ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
            
            ax.set_ylim([min_m_plot, max_m_plot])
            ax.grid(True, which='major', color='lightgray', linestyle='--', alpha=0.5)
            ax.grid(True, which='minor', color='whitesmoke', linestyle=':', alpha=0.3)
            ax.axhline(0.1, color='black', ls=':', alpha=0.8, lw=2.0)
            
            if t_cross_avg:
                t_med = np.median(t_cross_avg)
                ax.axvline(t_med, color='blue', ls='--', alpha=0.7, lw=2.0)
                
            if col == 0: ax.set_ylabel(r"Masa del embrión [$M_\oplus$]", fontsize=16)
            last_in_col = max([idx for idx in range(n_amp) if idx % cols == col])
            if i == last_in_col or i + cols >= len(axes): ax.set_xlabel("Tiempo [Myr]", fontsize=16)
            
        for j in range(n_amp, len(axes)):
            axes[j].set_visible(False)
            
        custom_lines = []
        labels = []
        for i, ng in enumerate(all_ngaps):
            c_val = norm(ng) if isinstance(ng, (int, float)) else norm(i)
            custom_lines.append(Line2D([0], [0], color=cmap(c_val), lw=3.5, ls=styles_list[i % len(styles_list)]))
            labels.append(f"{ng} Gaps" if isinstance(ng, (int, float)) else str(ng))
            
        custom_lines.append(Line2D([0], [0], color='black', lw=2.0, ls=':', alpha=0.8))
        labels.append("Umbral (0.1 $M_\oplus$)")
        if any('m_iso_e' in item for item in data):
            custom_lines.append(Line2D([0], [0], color='saddlebrown', lw=2.5, ls='--', alpha=0.7))
            labels.append("$M_{\mathrm{iso}}$")
            
        ncols = min(len(labels), 6)
        fig.legend(custom_lines, labels, loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=ncols, frameon=False, fontsize=14)
            
        fig.suptitle("Evolución Temporal por Amplitud (Mosaico)", fontsize=16, y=0.98)
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_amp_mosaic.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_lines_grouped_by_ngap_mosaic(data, fig_dir):
        import math
        if all('n_gaps' in d for d in data):
            group_key = 'n_gaps'
        else:
            group_key = 'gap_type'
            
        groups = {}
        for item in data:
            val = item.get(group_key)
            if val not in groups: groups[val] = []
            groups[val].append(item)
            
        ngaps = sorted(list(groups.keys()))
        if not ngaps: return
        n_ngap = len(ngaps)
        cols = min(n_ngap, 3)
        rows = math.ceil(n_ngap / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharex=True, sharey=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0, top=0.90)
        
        if not isinstance(axes, np.ndarray): axes = np.array([axes])
        axes = axes.flatten()
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list("custom_dist", ["darkred", "darkorange", "steelblue", "navy"])
        all_amps = sorted(list(set([item['amp_val'] for item in data])))
        from matplotlib.colors import Normalize
        norm = Normalize(vmin=min(all_amps), vmax=max(all_amps))
        
        min_m = np.nanmin([np.nanmin(item['mass_e']) for item in data])
        if min_m > 1e-2: min_m = 1e-2
        min_m_plot = min_m / 3.0
        max_m_plot = 10.0
        
        for i, (ax, ng) in enumerate(zip(axes, ngaps)):
            title_text = f"{ng} Gaps" if isinstance(ng, (int, float)) else str(ng)
            ax.text(0.05, 0.95, r"$N_{\rm gaps} = $ " + title_text, transform=ax.transAxes, 
                    ha='left', va='top', fontsize=14,
                    bbox=dict(facecolor='white', edgecolor='lightgray', alpha=0.8))
            
            ax.tick_params(direction='in', top=True, right=True, which='major', labelsize=14, length=6, width=1.5)
            ax.tick_params(direction='in', top=True, right=True, which='minor', length=3, width=0.5)
            
            runs_ng = sorted(groups[ng], key=lambda x: x['amp_val'])
            t_cross_avg = []
            
            styles_list = ['-', '--', '-.', ':', (0, (5, 5)), (0, (3, 5, 1, 5))]
            for item in runs_ng:
                color = cmap(norm(item['amp_val']))
                idx = all_amps.index(item['amp_val'])
                ls = styles_list[idx % len(styles_list)]
                
                t_myr = item['times_yr'] / 1e6
                mass_e = item['mass_e']
                if 'm_iso_e' in item:
                    ax.plot(t_myr, item['m_iso_e'], color='saddlebrown', ls='--', alpha=0.7, lw=2.0, zorder=1)
                
                mask = np.isfinite(t_myr) & np.isfinite(mass_e)
                t_valid = t_myr[mask]
                m_valid = mass_e[mask]
                if len(t_valid) > 0 and t_valid[-1] < 10.0:
                    t_valid = np.append(t_valid, 10.0)
                    m_valid = np.append(m_valid, m_valid[-1])
                    
                ax.plot(t_valid, m_valid, color=color, alpha=0.8, lw=3.5, zorder=2, ls=ls)
                if item['t_cross_1au']: t_cross_avg.append(item['t_cross_1au'] / 1e6)
                
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim([1e-3, 10.0])
            
            import matplotlib.ticker as ticker
            row, col = i // cols, i % cols
            if col == 0:
                ax.set_xticks([1e-3, 1e-2, 1e-1, 1e0, 1e1])
            else:
                ax.set_xticks([1e-2, 1e-1, 1e0, 1e1])
            ax.xaxis.set_major_formatter(ticker.LogFormatterMathtext())
            
            ax.set_ylim([min_m_plot, max_m_plot])
            ax.grid(True, which='major', color='lightgray', linestyle='--', alpha=0.5)
            ax.grid(True, which='minor', color='whitesmoke', linestyle=':', alpha=0.3)
            ax.axhline(0.1, color='black', ls=':', alpha=0.8, lw=2.0)
            
            if t_cross_avg:
                t_med = np.median(t_cross_avg)
                ax.axvline(t_med, color='blue', ls='--', alpha=0.7, lw=2.0)
                
            if col == 0: ax.set_ylabel(r"Masa del embrión [$M_\oplus$]", fontsize=16)
            last_in_col = max([idx for idx in range(n_ngap) if idx % cols == col])
            if i == last_in_col or i + cols >= len(axes): ax.set_xlabel("Tiempo [Myr]", fontsize=16)
            
        for j in range(n_ngap, len(axes)):
            axes[j].set_visible(False)
            
        custom_lines = [Line2D([0], [0], color=cmap(norm(amp)), lw=3.5, ls=styles_list[i % len(styles_list)]) for i, amp in enumerate(all_amps)]
        labels = [f"A = {amp}" for amp in all_amps]
        custom_lines.append(Line2D([0], [0], color='black', lw=2.0, ls=':', alpha=0.8))
        labels.append("Umbral (0.1 $M_\oplus$)")
        if any('m_iso_e' in item for item in data):
            custom_lines.append(Line2D([0], [0], color='saddlebrown', lw=2.5, ls='--', alpha=0.7))
            labels.append("$M_{\mathrm{iso}}$")
            
        ncols = min(len(labels), 6)
        fig.legend(custom_lines, labels, loc='lower center', bbox_to_anchor=(0.5, 0.92), ncol=ncols, frameon=False, fontsize=14)
            
        fig.suptitle(r"Evolución Temporal por $N_{\mathrm{gaps}}$ (Mosaico)", fontsize=16, y=0.98)
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_ngap_mosaic.png"), dpi=200, bbox_inches='tight')
        plt.close()

class PlotterSmooth:
    @staticmethod
    def plot_evolution(data, fig_dir, alpha_val=None):
        if not data: return
        
        fig, ax = plt.subplots(figsize=(10, 6))
        title = "Evolución Temporal - Disco Smooth"
        if alpha_val is not None: title += rf" ($\alpha = {alpha_val}$)"
        ax.set_title(title, fontsize=18)
        
        t_cross_avg = []
        
        for item in data:
            color = 'blue'
            completed = item.get('completed', True)
            ls = '-' if completed else '-.'
            
            t_myr_log = np.log10(item['times_yr'] / 1e6)
            mass_e_log = np.log10(item['mass_e'])
            
            if 'm_iso_e' in item:
                m_iso_log = np.log10(item['m_iso_e'])
                ax.plot(t_myr_log, m_iso_log, color=color, ls=':', alpha=0.4, lw=1.5, zorder=1)
            
            mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
            ax.plot(t_myr_log[mask], mass_e_log[mask], color=color, alpha=0.8, lw=3.0, ls=ls)
            
            if item['t_cross_1au']: t_cross_avg.append(np.log10(item['t_cross_1au'] / 1e6))
            
        ax.set_xlim([-3, 1])
        min_m_local = np.nanmin([np.nanmin(np.log10(item['mass_e'])) for item in data])
        max_m_local = np.nanmax([np.nanmax(np.log10(item['mass_e'])) for item in data])
        padding = (max_m_local - min_m_local) * 0.1
        if padding == 0: padding = 0.2
        ax.set_ylim([min_m_local - padding, max_m_local + padding])
        ax.grid(True, alpha=0.3)
        
        if t_cross_avg:
            t_med = np.median(t_cross_avg)
            ax.axvline(t_med, color='black', ls='--', alpha=0.6, lw=2.0, label="Snowline a 1au")
            
        ax.set_ylabel(r"Log(M_embrion [$M_\oplus$])")
        ax.set_xlabel("Log(t [Myr])")
        
        custom_lines = []
        labels = []
        
        if any(not item.get('completed', True) for item in data):
            custom_lines.append(Line2D([0], [0], color='blue', lw=3.0, ls='-.'))
            labels.append("Incompleta (< 10Myr)")
        else:
            custom_lines.append(Line2D([0], [0], color='blue', lw=3.0, ls='-'))
            labels.append("Evolución M_embr")
            
        if any('m_iso_e' in item for item in data):
            custom_lines.append(Line2D([0], [0], color='blue', lw=1.5, ls=':', alpha=0.4))
            labels.append("M_iso (límite)")
            
        if t_cross_avg:
            custom_lines.append(Line2D([0], [0], color='black', ls='--', alpha=0.6, lw=2.0))
            labels.append("Snowline a 1au")
            
        ax.legend(custom_lines, labels, loc="upper left", fontsize=11)
        
        plt.tight_layout()
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "evo_smooth.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_heatmaps(data, fig_dir):
        if not data: return
        
        alphas = sorted(list(set([d.get('alpha', 0.0) for d in data])))
        m0_vals = sorted(list(set([d.get('m0_earth', 0.0) for d in data])))
        
        if len(alphas) == 0 or len(m0_vals) == 0: return
        
        Z_mass = np.full((len(m0_vals), len(alphas)), np.nan)
        Z_h2o = np.full((len(m0_vals), len(alphas)), np.nan)
        
        for item in data:
            alpha = item.get('alpha', 0.0)
            m0 = item.get('m0_earth', 0.0)
            if alpha not in alphas or m0 not in m0_vals: continue
            
            i = m0_vals.index(m0)
            j = alphas.index(alpha)
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(alphas)), np.arange(len(m0_vals))
        
        width = max(9, len(alphas) * 1.3)
        height = max(7, len(m0_vals) * 1.1)
        
        fig, ax = plt.subplots(figsize=(width, height))
        vmax_mass = np.nanmax(Z_mass) if not np.all(np.isnan(Z_mass)) else 1.0
        ax.imshow(Z_mass, origin='lower', aspect='auto', cmap='viridis', norm=LogNorm(vmin=np.nanmin(Z_mass), vmax=vmax_mass))
        for i in range(len(m0_vals)):
            for j in range(len(alphas)):
                val = Z_mass[i, j]
                if not np.isnan(val):
                    color = 'white' if val < np.nanpercentile(Z_mass, 50) else 'black'
                    ax.text(j, i, rf"{val:.3f} $M_\oplus$", ha="center", va="center", color=color, fontsize=10)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(alphas)
        ax.set_yticks(y_idx)
        ax.set_yticklabels([f"{m0:g}" for m0 in m0_vals])
        ax.set_xlabel(r"Viscosidad ($\alpha$)")
        ax.set_ylabel(r"Masa Inicial del Embrión [$M_\oplus$]")
        ax.set_title("Masa Final del Embrión - Smooth")
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_masa.png"), dpi=200, bbox_inches='tight')
        plt.close()
        
        fig, ax = plt.subplots(figsize=(width, height))
        vmax_h2o = np.nanmax(Z_h2o) if not np.all(np.isnan(Z_h2o)) else 100.0
        if vmax_h2o == 0: vmax_h2o = 1.0
        ax.imshow(Z_h2o, origin='lower', aspect='auto', cmap='Blues', vmin=0, vmax=vmax_h2o)
        for i in range(len(m0_vals)):
            for j in range(len(alphas)):
                val = Z_h2o[i, j]
                if not np.isnan(val):
                    color = 'white' if val > (vmax_h2o * 0.6) else 'black'
                    ax.text(j, i, f"{val:.1f}%", ha="center", va="center", color=color, fontsize=11)
        ax.set_xticks(x_idx)
        ax.set_xticklabels(alphas)
        ax.set_yticks(y_idx)
        ax.set_yticklabels([f"{m0:g}" for m0 in m0_vals])
        ax.set_xlabel(r"Viscosidad ($\alpha$)")
        ax.set_ylabel(r"Masa Inicial del Embrión [$M_\oplus$]")
        ax.set_title("Fracción de Agua Final - Smooth")
        plt.savefig(os.path.join(fig_dir, "heatmap_h2o.png"), dpi=200, bbox_inches='tight')
        plt.close()

    @staticmethod
    def plot_categorical_map(data, fig_dir):
        if not data: return
        
        alphas = sorted(list(set([d.get('alpha', 0.0) for d in data])))
        m0_vals = sorted(list(set([d.get('m0_earth', 0.0) for d in data])))
        
        if len(alphas) == 0 or len(m0_vals) == 0: return
        
        Z_mass = np.full((len(m0_vals), len(alphas)), np.nan)
        Z_h2o = np.full((len(m0_vals), len(alphas)), np.nan)
        
        for item in data:
            alpha = item.get('alpha', 0.0)
            m0 = item.get('m0_earth', 0.0)
            if alpha not in alphas or m0 not in m0_vals: continue
            
            i = m0_vals.index(m0)
            j = alphas.index(alpha)
            Z_mass[i, j] = item['mass_e'][-1]
            Z_h2o[i, j] = item['frac_h2o_final']
            
        x_idx, y_idx = np.arange(len(alphas)), np.arange(len(m0_vals))
        
        width = max(10, len(alphas) * 1.3)
        height = max(8, len(m0_vals) * 1.1)
        fig, ax = plt.subplots(figsize=(width, height))
        
        ax.set_facecolor('#e0e0e0')
        Z_h2o_masked = np.ma.masked_where(Z_mass < 0.1, Z_h2o)
        
        cmap = plt.get_cmap('PuBu')
        vmax_h2o_universal = 50.0
        
        im = ax.imshow(Z_h2o_masked, origin='lower', aspect='auto', cmap=cmap, vmin=0, vmax=vmax_h2o_universal)
        
        cbar = fig.colorbar(im, ax=ax, pad=0.02)
        cbar.set_label(r"Fracción de Agua Final (%)", fontsize=14)
        
        for i in range(len(m0_vals)):
            for j in range(len(alphas)):
                mass = Z_mass[i, j]
                if not np.isnan(mass) and mass >= 0.1:
                    if i == 0 or np.isnan(Z_mass[i-1, j]) or Z_mass[i-1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i-0.5, i-0.5], color='black', lw=3)
                    if i == len(m0_vals)-1 or np.isnan(Z_mass[i+1, j]) or Z_mass[i+1, j] < 0.1:
                        ax.plot([j-0.5, j+0.5], [i+0.5, i+0.5], color='black', lw=3)
                    if j == 0 or np.isnan(Z_mass[i, j-1]) or Z_mass[i, j-1] < 0.1:
                        ax.plot([j-0.5, j-0.5], [i-0.5, i+0.5], color='black', lw=3)
                    if j == len(alphas)-1 or np.isnan(Z_mass[i, j+1]) or Z_mass[i, j+1] < 0.1:
                        ax.plot([j+0.5, j+0.5], [i-0.5, i+0.5], color='black', lw=3)
                        
        for i in range(len(m0_vals)):
            for j in range(len(alphas)):
                mass = Z_mass[i, j]
                water = Z_h2o[i, j]
                if not np.isnan(mass) and not np.isnan(water):
                    if mass >= 0.1:
                        text_color = 'white' if water > (vmax_h2o_universal * 0.5) else 'black'
                        label = f"{mass:.2f}\n{water:.1f}%"
                        ax.text(j, i, label, ha="center", va="center", color=text_color, fontsize=11, fontweight='bold')
                        
        ax.set_xticks(x_idx)
        ax.set_xticklabels(alphas)
        ax.set_yticks(y_idx)
        ax.set_yticklabels([f"{m0:g}" for m0 in m0_vals])
        ax.set_xlabel(r"Viscosidad ($\alpha$)", fontsize=14)
        ax.set_ylabel(r"Masa Inicial del Embrión [$M_\oplus$]", fontsize=14)
        ax.tick_params(labelsize=12)
        
        ax.annotate("Valores en celdas: Masa final [$M_\\oplus$] (arriba) | Fracción de Agua (abajo)", 
                    xy=(0.5, -0.15), xycoords='axes fraction', ha='center', va='top', fontsize=12)
        
        ax.set_title("Mapa de Evolución Smooth (Continente del Éxito)", fontsize=16)
        
        os.makedirs(fig_dir, exist_ok=True)
        plt.savefig(os.path.join(fig_dir, "heatmap_categorico.png"), dpi=200, bbox_inches='tight')
        plt.close()

import matplotlib.colors as colors

class PlotterBenchmarks:
    
    @staticmethod
    def _set_publication_style():
        """Aplica un tema de alta calidad tipo 'Astrophysical Journal' a las gráficas."""
        plt.rcParams.update({
            'font.family': 'serif',
            'font.size': 12,
            'axes.labelsize': 14,
            'axes.titlesize': 15,
            'axes.linewidth': 1.5,
            'xtick.labelsize': 12,
            'ytick.labelsize': 12,
            'xtick.major.size': 6,
            'xtick.minor.size': 3,
            'xtick.major.width': 1.5,
            'ytick.major.size': 6,
            'ytick.minor.size': 3,
            'ytick.major.width': 1.5,
            'xtick.direction': 'in',
            'ytick.direction': 'in',
            'legend.fontsize': 11,
            'legend.frameon': True,
            'legend.edgecolor': 'black',
            'figure.dpi': 300, # Alta resolución obligatoria
            'axes.grid': True,
            'grid.alpha': 0.3,
            'grid.linestyle': '--'
        })

    @staticmethod
    def plot_bubble_chart(df, fig_path, alpha_val="10^{-3}"):
        PlotterBenchmarks._set_publication_style()
        fig, ax = plt.subplots(figsize=(8, 6))
        
        # Mapa de color personalizado: Rojo (0%) -> Azul (10%) -> Morado (>=20%)
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.cm as cm
        custom_colors = ['#d73027', '#0077bb', '#660099'] # Rojo, Azul, Morado
        cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
        
        # 1. Separación Lógica de los Datos
        mask_success = df['M_final'] >= 0.1
        mask_fail = df['M_final'] < 0.1
        
        df_success = {k: v[mask_success] for k, v in df.items()}
        df_fail = {k: v[mask_fail] for k, v in df.items()}
        
        # 2. Graficar los fracasos (Triángulos grises pequeños)
        if len(df_fail['M_final']) > 0:
            ax.scatter(
                df_fail['r_gap'], df_fail['M_gap'],
                s=40, # Tamaño fijo y pequeño
                c='#d3d3d3', # Color gris claro
                marker='v', # Triángulo apuntando hacia abajo
                alpha=0.7,
                edgecolors='black',
                linewidth=0.5
            )
            
        # 3. Graficar los éxitos (Círculos escalados)
        scatter = None
        if len(df_success['M_final']) > 0:
            sizes_success = df_success['M_final'] * 400 + 20
            scatter = ax.scatter(
                df_success['r_gap'], df_success['M_gap'],
                s=sizes_success,
                c=df_success['frac_h2o'],
                cmap=cmap_custom,
                vmin=0.0, vmax=0.20,
                marker='o', # Círculo clásico
                alpha=0.85,
                edgecolors='black',
                linewidth=0.7
            )
        
        # Escala logarítmica para el Gap Depth
        ax.set_yscale('log')
        ax.set_ylim(0.005, 5.0)
        
        # Márgenes en X
        if len(df['r_gap']) > 0:
            ax.set_xlim(min(df['r_gap'])-2, max(df['r_gap'])+2)
        
        ax.set_xlabel("Posición del Gap [AU]")
        ax.set_ylabel(r"Profundidad del Gap [$M_{\rm Jup}$]")
        # Forzar visualización de números en el eje log
        ax.set_yticks([0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 5.0])
        from matplotlib.ticker import ScalarFormatter
        ax.yaxis.set_major_formatter(ScalarFormatter())
        
        # Colorbar (Manejando el caso donde no hay éxitos)
        if scatter is not None:
            cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        else:
            sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, pad=0.02)
            
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        # Texto indicativo de turbulencia y leyenda geométrica
        legend_text = (r"$\alpha = {}$".format(alpha_val) + "\n" +
                       r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        
        ax.text(0.05, 0.95, legend_text,
                transform=ax.transAxes, fontsize=11,
                verticalalignment='top',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        
        plt.tight_layout()
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)
    @staticmethod
    def plot_hovmoller(t_array, r_array, eps_matrix, r_ice_array, fig_path, run_name):
        PlotterBenchmarks._set_publication_style()
        fig, ax = plt.subplots(figsize=(10, 6))
        
        t_myr = t_array / (3.15576e7 * 1e6)
        
        # 'inferno' o 'magma' son el estándar en astronomía para mapas de densidad
        pcm = ax.pcolormesh(
            t_myr, r_array, eps_matrix.T,
            norm=colors.LogNorm(vmin=1e-4, vmax=1e-1),
            cmap='inferno', shading='nearest'
        )
        
        # Líneas más gruesas y con contorno sutil para que destaquen sobre fondos claros/oscuros
        ax.plot(t_myr, r_ice_array, color='cyan', linestyle='--', linewidth=2.0, 
                label=r'Snowline Dinámica $R_{\rm ice}(t)$')
        
        ax.axhline(1.0, color='white', linestyle=':', linewidth=1.5, alpha=0.8, label='Embrión (1 AU)')
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-3, 10.0)
        ax.set_ylim(min(r_array), max(r_array))
        
        ax.set_xlabel("Tiempo [Myr]")
        ax.set_ylabel("Distancia Radial [AU]")
        ax.set_title(rf"Evolución de la Trampa de Polvo ($\Sigma_d/\Sigma_g$) - {run_name}", pad=15)
        
        cbar = fig.colorbar(pcm, ax=ax, extend='both', pad=0.02)
        cbar.set_label(r"Relación Polvo-Gas $\epsilon$", rotation=270, labelpad=20)
        
        # Fondo oscuro para la leyenda para que no rompa el gráfico
        legend = ax.legend(loc="upper right", framealpha=0.8, facecolor='black', edgecolor='white')
        plt.setp(legend.get_texts(), color='w')
        
        plt.tight_layout()
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_pebble_flux(t_array, flux_1au, fig_path, run_name):
        PlotterBenchmarks._set_publication_style()
        fig, ax = plt.subplots(figsize=(8, 5))
        t_myr = t_array / (3.15576e7 * 1e6)
        
        # Usamos np.abs porque el flujo radial apunta hacia la estrella (es negativo)
        flux_m_earth_per_myr = np.abs(flux_1au) * (3.15576e7 * 1e6) / 5.9722e27
        
        # Línea principal gruesa
        ax.plot(t_myr, flux_m_earth_per_myr, color='#1f77b4', linewidth=2.5, zorder=3)
        
        # Sombrear el área bajo la curva le da "peso" al flujo másico
        ax.fill_between(t_myr, 1e-10, flux_m_earth_per_myr, color='#1f77b4', alpha=0.2, zorder=2)
        
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(1e-3, 10.0)
        
        # Ajuste dinámico de los límites Y para evitar inversiones o gráficos vacíos
        # Consideramos "flujo válido" todo lo mayor a 1e-6 M_earth/Myr (ruido numérico por debajo de esto)
        valid_flux = flux_m_earth_per_myr[flux_m_earth_per_myr > 1e-6]
        if len(valid_flux) > 0:
            y_min = np.min(valid_flux) * 0.5
            y_max = np.max(valid_flux) * 5.0
            ax.set_ylim(y_min, max(y_max, 1e-2))
        else:
            # Caso de "Asfixia" real, forzamos un gráfico vacío con límites sanos
            ax.set_ylim(1e-6, 1e-2)
        
        ax.set_xlabel("Tiempo [Myr]")
        ax.set_ylabel(r"Flujo Neto de Pebbles a 1 AU [$M_\oplus / {\rm Myr}$]")
        ax.set_title(f"Tasa de Inyección de Sólidos - {run_name}", pad=15)
        
        plt.tight_layout()
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_bubble_chart_mosaic(data_por_alpha, alphas_to_plot, fig_path):
        import math
        PlotterBenchmarks._set_publication_style()
        
        n_panels = len(alphas_to_plot)
        cols = min(n_panels, 3)
        rows = math.ceil(n_panels / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols + 1, 4.5 * rows), sharey=True, sharex=True)
        
        # 2. Espaciado Interno Cero
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.cm as cm
        custom_colors = ['#d73027', '#0077bb', '#660099'] 
        cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
        
        import numpy as np
        axes = np.atleast_1d(axes).flatten()
            
        for i in range(len(axes)):
            ax = axes[i]
            if i >= n_panels:
                ax.set_visible(False)
                continue
                
            alpha = alphas_to_plot[i]
            df = data_por_alpha.get(alpha, None)
            
            # Ticks Interiores y Espejados SIEMPRE
            ax.tick_params(direction='in', top=True, right=True, which='both')
            
            if df is None or len(df['M_final']) == 0:
                ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                        ha='center', va='top', fontsize=13, 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
                continue
                
            mask_success = df['M_final'] >= 0.1
            mask_fail = df['M_final'] < 0.1
            
            df_success = {k: v[mask_success] for k, v in df.items()}
            df_fail = {k: v[mask_fail] for k, v in df.items()}
            
            if len(df_fail['M_final']) > 0:
                ax.scatter(
                    df_fail['r_gap'], df_fail['M_gap'],
                    s=40, c='#d3d3d3', marker='v', alpha=0.7,
                    edgecolors='black', linewidth=0.5
                )
                
            if len(df_success['M_final']) > 0:
                sizes_success = df_success['M_final'] * 400 + 20
                ax.scatter(
                    df_success['r_gap'], df_success['M_gap'],
                    s=sizes_success, c=df_success['frac_h2o'],
                    cmap=cmap_custom, vmin=0.0, vmax=0.20,
                    marker='o', alpha=0.85,
                    edgecolors='black', linewidth=0.7
                )
            
            # El título lo ponemos adentro del gráfico para que no choque con hspace=0
            ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='center', va='top', fontsize=13, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
            
            # Solo la columna izquierda muestra Y label
            if i % cols == 0:
                ax.set_ylabel(r"Profundidad del Gap [$M_{\rm Jup}$]")
            
            # Solo la última fila muestra X label
            # Nota: si el panel está en la última fila (o si el panel abajo está oculto)
            if i >= cols * (rows - 1) or i + cols >= n_panels:
                ax.set_xlabel("Posición del Gap [AU]")
            
            # Escala logarítmica y límites
            ax.set_yscale('log')
            ax.set_ylim(0.005, 5.0)
            from matplotlib.ticker import ScalarFormatter
            ax.yaxis.set_major_formatter(ScalarFormatter())
            ax.set_yticks([0.01, 0.05, 0.1, 0.3, 0.5, 1.0, 3.0, 5.0])
            
            if len(df['r_gap']) > 0:
                ax.set_xlim(min(df['r_gap'])-2, max(df['r_gap'])+2)
        
        # 4. Elementos Globales (Colorbar y Leyendas)
        sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
        sm.set_array([])
        
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        # Leyenda maestra en el primer panel
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        axes[0].text(0.05, 0.05, legend_text,
                transform=axes[0].transAxes, fontsize=11,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)

    @staticmethod
    def plot_hovmoller_mosaic(data_hov, alphas_to_plot, fig_path, title=None, r_min_val=None):
        PlotterBenchmarks._set_publication_style()
        import math
        import matplotlib.colors as colors
        import numpy as np
        
        n_panels = len(alphas_to_plot)
        cols = min(n_panels, 3)
        rows = math.ceil(n_panels / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharey=True, sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        axes = np.atleast_1d(axes).flatten()
        
        pcm = None
        for i in range(len(axes)):
            ax = axes[i]
            if i >= n_panels:
                ax.set_visible(False)
                continue
                
            alpha = alphas_to_plot[i]
            if alpha not in data_hov:
                ax.set_title(rf"$\alpha = {alpha}$ (Sin datos)", pad=10)
                continue
                
            d = data_hov[alpha]
            t_myr = d['t_arr'] / (3.15576e7 * 1e6)
            
            pcm = ax.pcolormesh(
                t_myr, d['r_arr'], d['eps_mat'].T,
                norm=colors.LogNorm(vmin=1e-4, vmax=1e-1),
                cmap='inferno', shading='nearest'
            )
            
            r_ice_real = 2.73 * (t_myr / 0.2)**(-0.5)
            r_ice_real = np.clip(r_ice_real, 0.5, 2.73)
            
            ax.plot(t_myr, r_ice_real, color='cyan', linestyle='--', linewidth=2.0)
            ax.axhline(1.0, color='white', linestyle=':', linewidth=1.5, alpha=0.8)
            
            if r_min_val is not None:
                ax.axhline(r_min_val, color='lime', linestyle='-.', linewidth=1.5, alpha=0.8)
            
            ax.tick_params(direction='in', top=True, right=True, which='both')
            
            ax.text(0.05, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=13, color='white',
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5, edgecolor='none'))
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1e-3, 10.0)
            
            if r_min_val is not None:
                ax.set_ylim(r_min_val, 100.0)
            else:
                # Para el caso general/smooth, usualmente empieza en 0.5
                ax.set_ylim(0.5, 100.0)
                
            if i % cols == 0:
                ax.set_ylabel("Distancia Radial [AU]")
            if i >= cols * (rows - 1) or i + cols >= n_panels:
                ax.set_xlabel("Tiempo [Myr]")
                
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        if pcm is not None:
            cbar = fig.colorbar(pcm, cax=cbar_ax, extend='both')
            cbar.set_label(r"Relación Polvo-Gas $\epsilon$", rotation=270, labelpad=20)
        
        if title:
            fig.suptitle(title, fontsize=16, y=0.94)
            
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_pebble_flux_mosaic(data_hov, alphas_to_plot, fig_path):
        PlotterBenchmarks._set_publication_style()
        import math
        import numpy as np
        
        n_panels = len(alphas_to_plot)
        cols = min(n_panels, 3)
        rows = math.ceil(n_panels / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharey=True, sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        axes = np.atleast_1d(axes).flatten()
        
        for i in range(len(axes)):
            ax = axes[i]
            if i >= n_panels:
                ax.set_visible(False)
                continue
                
            alpha = alphas_to_plot[i]
            if alpha not in data_hov:
                continue
                
            d = data_hov[alpha]
            t_myr = d['t_arr'] / (3.15576e7 * 1e6)
            flux_m_earth_per_myr = np.abs(d['flux_arr']) * (3.15576e7 * 1e6) / 5.9722e27
            
            ax.plot(t_myr, flux_m_earth_per_myr, color='#1f77b4', linewidth=2.5, zorder=3)
            ax.fill_between(t_myr, 1e-10, flux_m_earth_per_myr, color='#1f77b4', alpha=0.2, zorder=2)
            
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlim(1e-3, 10.0)
            ax.set_ylim(1e-6, 1e4) # Ampliado para capturar las avalanchas de pebbles
            
            ax.tick_params(direction='in', top=True, right=True, which='both')
            
            ax.text(0.05, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=13, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
            
            if i % cols == 0:
                ax.set_ylabel(r"Flujo 1 AU [$M_\oplus / {\rm Myr}$]")
            if i >= cols * (rows - 1) or i + cols >= n_panels:
                ax.set_xlabel("Tiempo [Myr]")
                
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)

    @staticmethod
    def plot_flux_amax_mosaic(data_hov, alphas_to_plot, fig_path, title=None):
        PlotterBenchmarks._set_publication_style()
        import math
        import numpy as np
        
        n_panels = len(alphas_to_plot)
        cols = min(n_panels, 3)
        rows = math.ceil(n_panels / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 4 * rows), sharey=True, sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        axes = np.atleast_1d(axes).flatten()
        
        # Almacenar ejes gemelos para compartir la escala
        twin_axes = []
        for ax in axes:
            ax_twin = ax.twinx()
            twin_axes.append(ax_twin)
            
        for i in range(len(twin_axes)):
            if i > 0:
                twin_axes[i].sharey(twin_axes[0])
        
        for i in range(len(axes)):
            ax = axes[i]
            ax_twin = twin_axes[i]
            
            if i >= n_panels:
                ax.set_visible(False)
                ax_twin.set_visible(False)
                continue
                
            alpha = alphas_to_plot[i]
            if alpha not in data_hov:
                continue
                
            d = data_hov[alpha]
            t_myr = d['t_arr'] / (3.15576e7 * 1e6)
            flux_m_earth_per_myr = np.abs(d['flux_arr']) * (3.15576e7 * 1e6) / 5.9722e27
            amax_cm = d['amax_arr']
            
            # Máscara Física: Ocultar a_max cuando no hay polvo (flujo colapsado)
            # Threshold de 1e-4 M_earth/Myr es un límite inferior muy seguro.
            amax_cm_masked = np.where(flux_m_earth_per_myr > 1e-4, amax_cm, np.nan)
            
            # 1. Curva de Flujo (Eje Izquierdo)
            ax.plot(t_myr, flux_m_earth_per_myr, color='#1f77b4', linewidth=2.5, zorder=3, label="Flujo de Masa")
            ax.fill_between(t_myr, 1e-6, flux_m_earth_per_myr, color='#1f77b4', alpha=0.2, zorder=2)
            
            # 2. Curva de a_max (Eje Derecho Gemelo) usando la versión enmascarada
            ax_twin.plot(t_myr, amax_cm_masked, color='dimgray', linestyle='--', linewidth=2.0, zorder=4, label=r"$a_{\rm max}$")
            
            # Escalas y Límites
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax_twin.set_yscale('log')
            
            ax.set_xlim(1e-3, 10.0)
            ax.set_ylim(1e-6, 1e4)
            ax_twin.set_ylim(1e-4, 1e2)
            
            # Ticks hacia adentro
            ax.tick_params(direction='in', top=True, right=False, which='both')
            ax_twin.tick_params(direction='in', right=True, left=False, which='both')
            
            # Visibilidad de etiquetas en ejes compartidos
            if i % cols == 0:
                ax.set_ylabel(r"Flujo Neto a 1 AU [$M_\oplus/{\rm Myr}$]")
            else:
                ax.tick_params(labelleft=False)
                
            if i >= cols * (rows - 1) or i + cols >= n_panels:
                ax.set_xlabel("Tiempo [Myr]")
                
            if (i + 1) % cols == 0:
                ax_twin.set_ylabel(r"$a_{\rm max}$ [cm]", rotation=270, labelpad=20)
            else:
                ax_twin.tick_params(labelright=False)
            
            # Etiqueta del parámetro
            ax.text(0.05, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='left', va='top', fontsize=13, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'), zorder=5)
            
            # Leyenda unificada en el primer panel
            if i == 0:
                lines_1, labels_1 = ax.get_legend_handles_labels()
                lines_2, labels_2 = ax_twin.get_legend_handles_labels()
                ax.legend(lines_1 + lines_2, labels_1 + labels_2, loc='lower left', 
                          fontsize=11, framealpha=0.9, edgecolor='#cccccc')
                
        if title:
            fig.suptitle(title, fontsize=16, y=0.94)
            
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_smooth_global_benchmark(df, fig_path):
        PlotterBenchmarks._set_publication_style()
        fig, ax = plt.subplots(figsize=(8, 6))
        
        custom_colors = ['#d73027', '#0077bb', '#660099'] 
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.cm as cm
        cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
        
        v_frags = sorted(df['v_frag'].unique())
        alphas = sorted(df['alpha'].unique())
        
        v_map = {v: i for i, v in enumerate(v_frags)}
        a_map = {a: i for i, a in enumerate(alphas)}
        
        df_success = df[df['M_final'] >= 0.1]
        df_fail = df[df['M_final'] < 0.1]
        
        if len(df_fail) > 0:
            x_fail = [v_map[v] for v in df_fail['v_frag']]
            y_fail = [a_map[a] for a in df_fail['alpha']]
            ax.scatter(x_fail, y_fail, s=40, c='#d3d3d3', marker='v', alpha=0.7, edgecolors='black', linewidth=0.5)
            
        scatter = None
        if len(df_success) > 0:
            x_succ = [v_map[v] for v in df_success['v_frag']]
            y_succ = [a_map[a] for a in df_success['alpha']]
            sizes = df_success['M_final'] * 100 + 20
            frac_h2o = df_success['frac_h2o_percent'] / 100.0
            
            scatter = ax.scatter(x_succ, y_succ, s=sizes, c=frac_h2o, cmap=cmap_custom, vmin=0.0, vmax=0.20,
                                 marker='o', alpha=0.85, edgecolors='black', linewidth=0.7)
                                 
        ax.set_xticks(range(len(v_frags)))
        ax.set_xticklabels([f"{v}" for v in v_frags])
        ax.set_yticks(range(len(alphas)))
        ax.set_yticklabels([rf"${a}$" for a in alphas])
        
        ax.set_xlim(-0.5, len(v_frags) - 0.5)
        ax.set_ylim(-0.5, len(alphas) - 0.5)
        
        ax.set_xlabel(r"Velocidad de fragmentación $v_{\rm frag}$ [m/s]")
        ax.set_ylabel(r"Viscosidad $\alpha$")
        ax.set_title(r"Validación Global: Discos Smooth ($A=0$)", fontsize=15, pad=15)
        
        if scatter is not None:
            cbar = plt.colorbar(scatter, ax=ax, pad=0.02)
        else:
            sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
            sm.set_array([])
            cbar = plt.colorbar(sm, ax=ax, pad=0.02)
            
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        ax.text(0.05, 0.95, legend_text, transform=ax.transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
                
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_sinusoidal_benchmark(df, fig_path):
        PlotterBenchmarks._set_publication_style()
        
        v_frags = [1, 3, 10]
        alphas = sorted(df['alpha'].unique())
        amps = sorted(df['amp_val'].unique())
        
        a_map = {a: i for i, a in enumerate(alphas)}
        amp_map = {am: i for i, am in enumerate(amps)}
        
        fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)
        plt.subplots_adjust(wspace=0.0)
        
        custom_colors = ['#d73027', '#0077bb', '#660099'] 
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.cm as cm
        cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
        
        for i, vf in enumerate(v_frags):
            ax = axes[i]
            df_vf = df[df['v_frag'] == vf]
            
            ax.set_title(rf"$v_{{\rm frag}} = {vf}$ m/s", fontsize=14)
            ax.tick_params(direction='in', top=True, right=True, which='both')
            
            ax.set_xticks(range(len(amps)))
            ax.set_xticklabels([f"{am}" for am in amps])
            
            if i == 0:
                ax.set_yticks(range(len(alphas)))
                ax.set_yticklabels([rf"${a}$" for a in alphas])
                ax.set_ylabel(r"Viscosidad $\alpha$")
            
            ax.set_xlabel(r"Amplitud $A$")
            
            ax.set_xlim(-0.5, len(amps) - 0.5)
            ax.set_ylim(-0.5, len(alphas) - 0.5)
            
            if len(df_vf) == 0:
                continue
                
            df_success = df_vf[df_vf['M_final'] >= 0.1]
            df_fail = df_vf[df_vf['M_final'] < 0.1]
            
            if len(df_fail) > 0:
                x_fail = [amp_map[am] for am in df_fail['amp_val']]
                y_fail = [a_map[a] for a in df_fail['alpha']]
                ax.scatter(x_fail, y_fail, s=40, c='#d3d3d3', marker='v', alpha=0.7, edgecolors='black', linewidth=0.5)
                
            if len(df_success) > 0:
                x_succ = [amp_map[am] for am in df_success['amp_val']]
                y_succ = [a_map[a] for a in df_success['alpha']]
                sizes = df_success['M_final'] * 100 + 20
                frac_h2o = df_success['frac_h2o_percent'] / 100.0
                
                ax.scatter(x_succ, y_succ, s=sizes, c=frac_h2o, cmap=cmap_custom, vmin=0.0, vmax=0.20,
                           marker='o', alpha=0.85, edgecolors='black', linewidth=0.7)
                           
        fig.subplots_adjust(right=0.90)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.7])
        sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
        sm.set_array([])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        axes[0].text(0.05, 0.95, legend_text, transform=axes[0].transAxes, fontsize=11,
                verticalalignment='top', bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        
        fig.suptitle(r"Evolución en Discos Sinusoidales Intermedios ($N_{\rm gaps} = 5$)", fontsize=16, y=1.02)
        
        plt.savefig(fig_path, dpi=300, bbox_inches='tight')
        plt.close(fig)

    @staticmethod
    def plot_sinusoidal_bubble_mosaic(data_por_alpha, alphas_to_plot, fig_path):
        import math
        PlotterBenchmarks._set_publication_style()
        
        n_panels = len(alphas_to_plot)
        cols = min(n_panels, 3)
        rows = math.ceil(n_panels / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(4 * cols + 1, 4.5 * rows), sharey=True, sharex=True)
        plt.subplots_adjust(wspace=0.0, hspace=0.0)
        
        from matplotlib.colors import LinearSegmentedColormap
        import matplotlib.cm as cm
        custom_colors = ['#d73027', '#0077bb', '#660099'] 
        cmap_custom = LinearSegmentedColormap.from_list('water_fraction', custom_colors)
        
        import numpy as np
        axes = np.atleast_1d(axes).flatten()
            
        for i in range(len(axes)):
            ax = axes[i]
            if i >= n_panels:
                ax.set_visible(False)
                continue
                
            alpha = alphas_to_plot[i]
            df = data_por_alpha.get(alpha, None)
            
            ax.tick_params(direction='in', top=True, right=True, which='both')
            
            if df is None or len(df['M_final']) == 0:
                ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                        ha='center', va='top', fontsize=13, 
                        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
                continue
                
            mask_success = df['M_final'] >= 0.1
            mask_fail = df['M_final'] < 0.1
            
            df_success = {k: v[mask_success] for k, v in df.items()}
            df_fail = {k: v[mask_fail] for k, v in df.items()}
            
            if len(df_fail['M_final']) > 0:
                ax.scatter(
                    df_fail['amp_val'], df_fail['n_gaps'],
                    s=40, c='#d3d3d3', marker='v', alpha=0.7,
                    edgecolors='black', linewidth=0.5
                )
                
            if len(df_success['M_final']) > 0:
                sizes_success = df_success['M_final'] * 400 + 20
                ax.scatter(
                    df_success['amp_val'], df_success['n_gaps'],
                    s=sizes_success, c=df_success['frac_h2o'],
                    cmap=cmap_custom, vmin=0.0, vmax=0.20,
                    marker='o', alpha=0.85,
                    edgecolors='black', linewidth=0.7
                )
            
            ax.text(0.5, 0.95, rf"$\alpha = {alpha}$", transform=ax.transAxes, 
                    ha='center', va='top', fontsize=13, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='none'))
            
            if i % cols == 0:
                ax.set_ylabel(r"Número de Gaps ($N_{\rm gaps}$)")
            
            if i >= cols * (rows - 1) or i + cols >= n_panels:
                ax.set_xlabel(r"Amplitud ($A$)")
            
            ax.set_xlim(-0.5, 6.0)
            ax.set_ylim(-1, 26)
            
        sm = cm.ScalarMappable(cmap=cmap_custom, norm=plt.Normalize(vmin=0.0, vmax=0.20))
        sm.set_array([])
        
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        cbar = fig.colorbar(sm, cax=cbar_ax)
        cbar.set_label("Fracción de Agua Final")
        cbar.set_ticks([0.0, 0.1, 0.2])
        cbar.set_ticklabels(['0%', '10%', r'$\geq 20\%$'])
        
        legend_text = (r"$\bigcirc \geq 0.1 M_\oplus$ (Éxito)" + "\n" +
                       r"$\nabla < 0.1 M_\oplus$ (Fracaso)")
        axes[0].text(0.05, 0.05, legend_text,
                transform=axes[0].transAxes, fontsize=11,
                verticalalignment='bottom',
                bbox=dict(boxstyle='round,pad=0.4', facecolor='white', alpha=0.9, edgecolor='#cccccc'))
        
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        fig.savefig(fig_path, dpi=300)
        plt.close(fig)


class PlotterPopulation:
    @staticmethod
    def _set_publication_style():
        import matplotlib.pyplot as plt
        plt.style.use('default')
        plt.rcParams.update({
            'font.family': 'serif',
            'font.size': 16,
            'axes.labelsize': 20,
            'axes.titlesize': 22,
            'xtick.labelsize': 16,
            'ytick.labelsize': 16,
            'legend.fontsize': 14,
            'figure.dpi': 300,
            'savefig.dpi': 300,
            'xtick.direction': 'in',
            'ytick.direction': 'in',
            'xtick.top': True,
            'ytick.right': True,
            'xtick.minor.visible': True,
            'ytick.minor.visible': True,
        })

    @staticmethod
    def plot_poblacion_sintetica(data_list, fig_path, scenario_name=None, mass_threshold=1.0, water_threshold=0.10, color_metric="alpha"):
        import pandas as pd
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import matplotlib.ticker as ticker
        import matplotlib.patheffects as pe
        import numpy as np
        import os

        # Make sure directory exists
        os.makedirs(os.path.dirname(fig_path), exist_ok=True)

        PlotterPopulation._set_publication_style()

        rows = []
        for d in data_list:
            m_final = None
            f_water_final = None
            
            # Extract final mass handling multiple formats
            if 'm_emb' in d: m_final = d['m_emb'][-1]
            elif 'mass_e' in d: m_final = d['mass_e'][-1]
            elif 'M_final' in d: m_final = d['M_final']
            
            # Extract final water handling multiple formats
            if 'f_water' in d: f_water_final = d['f_water'][-1]
            elif 'frac_h2o_final' in d: f_water_final = d['frac_h2o_final'] / 100.0
            elif 'frac_h2o_percent' in d: f_water_final = d['frac_h2o_percent'] / 100.0
            
            if m_final is None or f_water_final is None:
                continue
                
            alpha = d.get('alpha', np.nan)
            vfrag = d.get('vfrag', np.nan)
            m0 = d.get('m0_earth', np.nan)

            f_post = 0.0
            if color_metric == 'f_post':
                t = np.array(d.get('times_yr', []))
                m = np.array(d.get('mass_e', []))
                t_c = d.get('t_cross_1au', 0)
                if t_c is None:
                    t_c = 0
                if len(t) > 0 and len(m) > 0:
                    m_c = m[t <= t_c][-1] if any(t <= t_c) else m[0]
                    f_post = (m[-1] - m_c) / m[-1] if m[-1] > 0 else 0

            rows.append({
                "alpha": alpha,
                "vfrag": vfrag,
                "m_embr0": m0,
                "f_water_final": f_water_final,
                "m_final": m_final,
                "f_post": f_post
            })

        df = pd.DataFrame(rows)
        if df.empty:
            print("No hay datos para graficar la población sintética.")
            return None, None

        # Filtrar masas iniciales para solo incluir desde 1e-1 hasta 1e-4
        valid_m0 = [m for m in df["m_embr0"].dropna().unique() if m >= 0.0001]

        # Mosaico 3x4: Filas = vfrag, Columnas = m_embr0
        cols_panels = sorted(valid_m0) # m_embr0 (de menor a mayor)
        rows_panels = sorted(df["vfrag"].dropna().unique()) # vfrag
        
        ncols = len(cols_panels)
        nrows = len(rows_panels)
        
        if ncols == 0 or nrows == 0:
            print("No hay valores validos para graficar.")
            return None, None

        # Tamaño de figura ajustado para 3x4 (ancho x alto)
        fig, axes = plt.subplots(nrows, ncols, figsize=(5.5 * ncols, 5.0 * nrows), sharex=True, sharey=True)
        
        if nrows == 1 and ncols == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = np.array([axes])
        elif ncols == 1:
            axes = np.array([[ax] for ax in axes])
            
        # Un espaciado mínimo en vertical evita que los números del eje Y (ej. 10^-3 y 10^1) choquen
        plt.subplots_adjust(wspace=0.0, hspace=0.08, right=0.85, top=0.95)

        if color_metric == "alpha":
            # Mapeo categórico de Alpha (Opción 1 del revisor)
            unique_alphas = sorted(df["alpha"].dropna().unique())
            alpha_to_idx = {a: i for i, a in enumerate(unique_alphas)}
            cmap = plt.get_cmap("turbo", len(unique_alphas))
            vmin, vmax = -0.5, len(unique_alphas) - 0.5
        else:
            # f_post continuo
            cmap = plt.get_cmap("magma")
            vmin, vmax = 0.0, 1.0

        sc = None
        for i, pval in enumerate(rows_panels):  # iterar vfrag en filas
            for j, mval in enumerate(cols_panels):  # iterar M0 en columnas
                ax = axes[i, j]
                sub = df[(df["vfrag"] == pval) & (df["m_embr0"] == mval)]
                
                if not sub.empty:
                    if color_metric == "alpha":
                        c_vals = sub["alpha"].map(alpha_to_idx)
                    else:
                        c_vals = sub["f_post"]
                        
                    sc = ax.scatter(
                        sub["f_water_final"], sub["m_final"],
                        c=c_vals, cmap=cmap, vmin=vmin, vmax=vmax,
                        s=80, alpha=0.7, edgecolors="none", zorder=3
                    )
                
                ax.set_yscale("log")
                ax.minorticks_on()
                ax.yaxis.set_minor_locator(ticker.LogLocator(base=10.0, subs=np.arange(2, 10) * 1.0, numticks=15))
                ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
                
                # Forzar visibilidad y tamaño de ticks mayores y menores
                ax.tick_params(which='major', direction='in', top=True, right=True, bottom=True, left=True, length=8, width=1.0)
                ax.tick_params(which='minor', direction='in', top=True, right=True, bottom=True, left=True, length=4, width=0.8)
                
                # Limite X dinámico dependiendo del max agua, minimo 0.3
                max_w = df["f_water_final"].max()
                x_limit = max(0.3, min(1.05, max_w + 0.1))
                ax.set_xlim(-0.01, x_limit)
                
                y_top = max(10.0, df["m_final"].max() * 2)
                ax.set_ylim(1e-3, y_top)

                # --- Sombreado de regiones de habitabilidad (Poblaciones Sintéticas) ---
                # Earth-like (0.05% a 10% agua, M >= mass_threshold)
                ax.fill_between([0.0005, water_threshold], mass_threshold, y_top, color='mediumseagreen', alpha=0.15, zorder=1)
                
                # Waterworlds (>10% agua, M >= mass_threshold)
                ax.fill_between([water_threshold, x_limit + 0.5], mass_threshold, y_top, color='dodgerblue', alpha=0.15, zorder=1)

                # Cuadrante umbral
                ax.axhline(mass_threshold, color="dimgray", linestyle="--", linewidth=1.5, alpha=0.5, zorder=2)
                ax.axvline(water_threshold, color="dimgray", linestyle="--", linewidth=1.5, alpha=0.5, zorder=2)

                # Etiquetas sutiles de las líneas y zonas (solo en el panel superior izquierdo)
                if i == 0 and j == 0:
                    text_path_effect = [pe.withStroke(linewidth=3, foreground="white")]
                    # Líneas guía
                    ax.text(0.01, mass_threshold * 0.8, rf"${mass_threshold:g} M_\oplus$ (Límite Supervivencia)", color="black", fontsize=12, alpha=0.9, ha="left", va="top", path_effects=text_path_effect)
                    ax.text(water_threshold + 0.01, 1e-3 * 1.5, f"{int(water_threshold*100)}% $H_2O$", color="black", fontsize=12, alpha=0.9, ha="left", va="bottom", rotation=90, path_effects=text_path_effect)
                    
                    # Nombres de Zonas
                    ax.text((0.0005 + water_threshold) / 2.0, y_top * 0.25, "Análogos\nTerrestres\n(0.05% - 10%)", color="black", fontsize=14, alpha=0.9, ha="center", va="center", fontweight="bold", path_effects=text_path_effect)
                    ax.text((water_threshold + x_limit) / 2.0, y_top * 0.25, "Mundos\nAcuáticos\n(> 10%)", color="black", fontsize=14, alpha=0.9, ha="center", va="center", fontweight="bold", path_effects=text_path_effect)

                # Títulos en la primera fila (M0)
                if i == 0:
                    exp_val = int(np.round(np.log10(mval)))
                    ax.set_title(rf"$M_{{emb,0}} = 10^{{{exp_val}}} M_\oplus$", pad=15)
                
                # Etiquetas X en la última fila
                if i == nrows - 1:
                    ax.set_xlabel("Fracción de Agua Final")
                    
                # Etiquetas de fila (vfrag) en la última columna
                if j == ncols - 1:
                    ax2 = ax.twinx()
                    ax2.set_yticks([])
                    ax2.set_ylabel(rf"$v_{{frag}} = {int(pval)}$ m/s", rotation=270, labelpad=30, fontsize=18)

        # Etiqueta Y centrada
        for i in range(nrows):
            if i == nrows // 2:
                axes[i, 0].set_ylabel(r"Masa Final del Planeta ($M_\oplus$)")

        if sc is not None:
            cbar_ax = fig.add_axes([0.88, 0.15, 0.035, 0.7])  # Más ancha
            if color_metric == "alpha":
                # Ticks en los enteros para centrar las etiquetas en los bloques de color
                cbar = fig.colorbar(sc, cax=cbar_ax, ticks=range(len(unique_alphas)))
                cbar.set_label(r"Viscosidad ($\alpha$)", fontsize=20)
                # Reemplazar los números enteros con los verdaderos valores de alpha
                cbar.ax.set_yticklabels([f"{a:g}" for a in unique_alphas])
                cbar.ax.tick_params(which='major', labelsize=16, length=0) # Ocultar las marcas físicas para look categórico
                cbar.ax.tick_params(which='minor', left=False, right=False) # Quitar minors porque es categórica
            else:
                cbar = fig.colorbar(sc, cax=cbar_ax)
                cbar.set_label(r"Masa post-hielo ($f_{\rm post}$)", fontsize=20)
                cbar.ax.tick_params(which='major', labelsize=16, length=8, direction='in')
                cbar.ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
                cbar.ax.tick_params(which='minor', length=4, direction='in')

        if fig_path:
            fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"-> Money Plot guardado en: {fig_path}")
        return fig, axes

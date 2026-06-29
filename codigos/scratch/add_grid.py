with open('C:/astro/Codigos practica + docs + papers/codigos/pipeline_analysis/plotters.py', 'r', encoding='utf-8') as f:
    content = f.read()

new_method = """
    @staticmethod
    def plot_grilla_cuantitativa(data_list, fig_path, scenario_name=None):
        import pandas as pd
        import numpy as np
        import matplotlib.pyplot as plt
        import seaborn as sns
        import os

        os.makedirs(os.path.dirname(fig_path), exist_ok=True)
        PlotterPopulation._set_publication_style()

        rows = []
        for d in data_list:
            m_final = None
            f_water_final = None
            if 'm_emb' in d: m_final = d['m_emb'][-1]
            elif 'mass_e' in d: m_final = d['mass_e'][-1]
            elif 'M_final' in d: m_final = d['M_final']
            
            if 'f_water' in d: f_water_final = d['f_water'][-1]
            elif 'frac_h2o_final' in d: f_water_final = d['frac_h2o_final'] / 100.0
            elif 'frac_h2o_percent' in d: f_water_final = d['frac_h2o_percent'] / 100.0
            
            m_iso = d.get('m_iso_e', 1000.0)
            if isinstance(m_iso, (list, np.ndarray)):
                m_iso = m_iso[-1] if len(m_iso) > 0 else 1000.0
                
            if m_final is None or f_water_final is None: continue
            
            rows.append({
                "alpha": d.get('alpha', np.nan),
                "v_frag": d.get('vfrag', np.nan),
                "M_c": m_final,
                "f_water": f_water_final,
                "M_iso": m_iso
            })
            
        df = pd.DataFrame(rows)
        if df.empty: return None

        df_cat = df.copy()
        
        # CATEGORIAS DE MASA
        df_cat['Cat_Masa'] = '< 0.1 M_E'
        mask_intermedio = (df_cat['M_c'] >= 0.1) & (df_cat['M_c'] < (df_cat['M_iso'] * 0.99))
        mask_aislado = (df_cat['M_c'] >= (df_cat['M_iso'] * 0.99))
        df_cat.loc[mask_intermedio, 'Cat_Masa'] = '> 0.1 M_E (No iso)'
        df_cat.loc[mask_aislado, 'Cat_Masa'] = 'Isolation Mass'

        # CATEGORIAS DE AGUA
        df_cat['Cat_Agua'] = '0% Agua'
        mask_baja = (df_cat['f_water'] > 0.0) & (df_cat['f_water'] <= 0.001)
        mask_trans = (df_cat['f_water'] > 0.001) & (df_cat['f_water'] <= 0.1)
        mask_ww = (df_cat['f_water'] > 0.1)
        
        df_cat.loc[mask_baja, 'Cat_Agua'] = 'Baja Agua (>0-0.1%)'
        df_cat.loc[mask_trans, 'Cat_Agua'] = 'Transición (0.1-10%)'
        df_cat.loc[mask_ww, 'Cat_Agua'] = 'Water Worlds (>10%)'
        
        df_cat['Exitoso'] = df_cat['M_c'] >= 0.1
        
        alphas = sorted(df_cat['alpha'].dropna().unique())
        vfrags = sorted(df_cat['v_frag'].dropna().unique(), reverse=True)
        
        matriz_exito = np.zeros((len(vfrags), len(alphas)))
        matriz_textos = [["" for _ in range(len(alphas))] for _ in range(len(vfrags))]
        
        for i, vfrag in enumerate(vfrags):
            for j, alpha in enumerate(alphas):
                subset = df_cat[(df_cat['alpha'] == alpha) & (df_cat['v_frag'] == vfrag)]
                total = len(subset)
                
                if total == 0:
                    continue
                    
                pct_exito = (subset['Exitoso'].sum() / total) * 100
                pct_iso = (len(subset[subset['Cat_Masa'] == 'Isolation Mass']) / total) * 100
                
                subset_exitosos = subset[subset['Exitoso']]
                total_exitosos = len(subset_exitosos) if len(subset_exitosos) > 0 else 1
                
                pct_secos = (len(subset_exitosos[subset_exitosos['Cat_Agua'] == '0% Agua']) / total_exitosos) * 100
                pct_baja = (len(subset_exitosos[subset_exitosos['Cat_Agua'] == 'Baja Agua (>0-0.1%)']) / total_exitosos) * 100
                pct_ww = (len(subset_exitosos[subset_exitosos['Cat_Agua'] == 'Water Worlds (>10%)']) / total_exitosos) * 100
                
                matriz_exito[i, j] = pct_exito
                
                texto = (
                    f"Llegan a >0.1 M_E: {pct_exito:.1f}%\\n"
                    f"Llegan a M_iso: {pct_iso:.1f}%\\n"
                    f"----\\n"
                    f"De los exitosos:\\n"
                    f"Secos (0%): {pct_secos:.1f}%\\n"
                    f"Baja Agua: {pct_baja:.1f}%\\n"
                    f"Waterworlds: {pct_ww:.1f}%"
                )
                matriz_textos[i][j] = texto

        fig, ax = plt.subplots(figsize=(14, 8))
        
        sns.heatmap(matriz_exito, annot=np.array(matriz_textos), fmt="", cmap="Blues", 
                    cbar_kws={'label': '% de planetas que superan 0.1 $M_\oplus$'}, 
                    xticklabels=[f"{a:g}" for a in alphas], yticklabels=[f"{v:g}" for v in vfrags], ax=ax, 
                    annot_kws={"size": 10, "family": "monospace"})
        
        ax.set_xlabel(r'Viscosidad Turbulenta ($\\alpha$)', fontsize=14)
        ax.set_ylabel(r'Velocidad de Fragmentación ($v_{frag}$) [m/s]', fontsize=14)
        title = 'Síntesis Cuantitativa de Poblaciones Planetarias'
        if scenario_name: title += f" - {scenario_name}"
        ax.set_title(title, fontsize=16)
        
        fig.tight_layout()
        fig.savefig(fig_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"-> Grilla Cuantitativa guardada en: {fig_path}")
        return fig, ax
"""

content = content + "\n" + new_method + "\n"

with open('C:/astro/Codigos practica + docs + papers/codigos/pipeline_analysis/plotters.py', 'w', encoding='utf-8') as f:
    f.write(content)

print("Appended grid method to plotters.py")

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.preprocessing import StandardScaler

# Set publication-style aesthetics
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "axes.labelsize": 12,
    "font.size": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
})

def plot_scatter_binned(df, ax, title):
    # Scatter transparente
    sns.scatterplot(
        data=df, x='log10_Memb0', y='frac_h2o_percent',
        hue='alpha', palette='viridis', alpha=0.3, ax=ax,
        edgecolor=None, s=30, legend=False
    )
    
    # Binning y medianas (agrupando por Memb0 exacto ya que son discretos: 1e-5, 1e-4...)
    # Agrupamos por alpha y log10_Memb0 para calcular medianas
    medians = df.groupby(['alpha', 'log10_Memb0'])['frac_h2o_percent'].median().reset_index()
    
    # Graficar las líneas de mediana
    sns.lineplot(
        data=medians, x='log10_Memb0', y='frac_h2o_percent',
        hue='alpha', palette='viridis', ax=ax, marker='o', 
        linewidth=2, markersize=8, legend='full'
    )
    ax.set_title(title)
    ax.set_xlabel(r"$\log_{10}(M_{\rm emb, 0}/M_\oplus)$")
    ax.set_ylabel(r"Fracción de Agua Final (%)")
    
    # Simplificar leyenda
    handles, labels = ax.get_legend_handles_labels()
    # Tomamos solo la mitad correspondiente al lineplot para no duplicar
    if len(handles) > 0:
        ax.legend(handles, labels, title=r"$\alpha$", loc='upper right', bbox_to_anchor=(1.0, 1.0))

def main():
    print("Cargando datos...")
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    
    print("Pre-procesando variables...")
    # Calcular variables logarítmicas
    # Asegurarnos de que M_emb0 y alpha > 0
    df = df[(df['M_emb0'] > 0) & (df['alpha'] > 0)].copy()
    
    df['log10_Memb0'] = np.log10(df['M_emb0'])
    df['log10_alpha'] = np.log10(df['alpha'])
    
    os.makedirs("data/benchmarks", exist_ok=True)
    report_path = "data/estadisticas_agua_reporte.txt"
    with open(report_path, "w", encoding='utf-8') as rep:
        rep.write("==================================================\n")
        rep.write("   REPORTE DE INFERENCIA ESTADÍSTICA (WATERWORLD)\n")
        rep.write("==================================================\n\n")
        
        rep.write(f"Total de simulaciones procesadas: {len(df)}\n\n")
        
        # 1. Correlaciones Simples
        print("Calculando correlaciones...")
        p_corr, p_pval = pearsonr(df['log10_Memb0'], df['frac_h2o_percent'])
        s_corr, s_pval = spearmanr(df['log10_Memb0'], df['frac_h2o_percent'])
        
        rep.write("--- 1. CORRELACIONES BIVARIADAS ---\n")
        rep.write(f"Correlación Pearson (Lineal):  r = {p_corr:.4f} (p-value = {p_pval:.2e})\n")
        rep.write(f"Correlación Spearman (Monótona): rho = {s_corr:.4f} (p-value = {s_pval:.2e})\n\n")

        # ====================================================
        # FIGURA 1: SCATTER GLOBAL
        # ====================================================
        print("Generando Figura 1: Scatter Global...")
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_scatter_binned(df, ax, "Global: Todas las $v_{\\rm frag}$ y Gaps")
        fig.tight_layout()
        fig.savefig("data/benchmarks/scatter_global_agua.png", dpi=300)
        plt.close(fig)

        # ====================================================
        # FIGURA 2: PANELES POR V_FRAG
        # ====================================================
        print("Generando Figura 2: Paneles por v_frag...")
        unique_v = sorted(df['v_frag'].unique())
        fig, axes = plt.subplots(1, len(unique_v), figsize=(6 * len(unique_v), 5), sharey=True)
        if len(unique_v) == 1:
            axes = [axes]
            
        for ax, v in zip(axes, unique_v):
            df_v = df[df['v_frag'] == v]
            plot_scatter_binned(df_v, ax, f"$v_{{\\rm frag}} = {v}$ m/s")
            
        fig.tight_layout()
        fig.savefig("data/benchmarks/scatter_vfrag_agua.png", dpi=300)
        plt.close(fig)

        # ====================================================
        # REGRESIONES OLS MÚLTIPLES
        # ====================================================
        print("Calculando modelos OLS estandarizados...")
        # Estandarizar variables independientes para comparar coeficientes
        scaler = StandardScaler()
        vars_to_scale = ['log10_Memb0', 'log10_alpha', 'v_frag', 'M_gap', 'r_gap']
        df_scaled = df.copy()
        df_scaled[vars_to_scale] = scaler.fit_transform(df[vars_to_scale])
        # No estandarizamos la dependiente para mantener la interpretabilidad en % de agua
        
        # Modelo 1: Simple Lineal
        model_simple = smf.ols("frac_h2o_percent ~ log10_Memb0 + log10_alpha + v_frag + M_gap + r_gap", data=df_scaled).fit()
        
        # Modelo 2: Con Interacción
        model_ext = smf.ols("frac_h2o_percent ~ log10_Memb0 * log10_alpha + v_frag + M_gap + r_gap", data=df_scaled).fit()
        
        rep.write("--- 2. REGRESIÓN MÚLTIPLE OLS (Variables Estandarizadas) ---\n")
        rep.write("Nota: Los coeficientes indican el cambio en % de agua por cada 1 desviación estándar de la variable.\n\n")
        
        rep.write("[MODELO SIMPLE]\n")
        rep.write(f"R-squared: {model_simple.rsquared:.4f}\n")
        rep.write(f"AIC: {model_simple.aic:.1f} | BIC: {model_simple.bic:.1f}\n")
        rep.write(model_simple.summary().tables[1].as_text())
        rep.write("\n\n")
        
        rep.write("[MODELO EXTENDIDO (Con Interacción)]\n")
        rep.write(f"R-squared: {model_ext.rsquared:.4f}\n")
        rep.write(f"AIC: {model_ext.aic:.1f} | BIC: {model_ext.bic:.1f}\n")
        rep.write(model_ext.summary().tables[1].as_text())
        rep.write("\n\n")
        
        if model_ext.aic < model_simple.aic:
            rep.write(f"-> CONCLUSIÓN AIC: El Modelo Extendido es superior (AIC menor por {model_simple.aic - model_ext.aic:.1f} puntos).\n\n")
        else:
            rep.write(f"-> CONCLUSIÓN AIC: El Modelo Simple es superior (AIC menor por {model_ext.aic - model_simple.aic:.1f} puntos).\n\n")

        # ====================================================
        # MACHINE LEARNING: RANDOM FOREST IMPORTANCE
        # ====================================================
        print("Entrenando Random Forest para Permutation Importance...")
        X = df_scaled[vars_to_scale]
        y = df['frac_h2o_percent']
        
        rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X, y)
        
        # Permutation Importance
        result = permutation_importance(rf, X, y, n_repeats=10, random_state=42, n_jobs=-1)
        importances = result.importances_mean
        
        rep.write("--- 3. IMPORTANCIA RELATIVA (Random Forest Permutation Importance) ---\n")
        rep.write("Mide la pérdida de precisión del modelo cuando la variable se mezcla al azar. A mayor valor, más poder predictivo.\n\n")
        
        imp_df = pd.DataFrame({'Variable': vars_to_scale, 'Importancia': importances})
        imp_df = imp_df.sort_values(by='Importancia', ascending=False)
        
        for idx, row in imp_df.iterrows():
            rep.write(f"{row['Variable']:>15}: {row['Importancia']:.4f}\n")
            
        rep.write("\n==================================================\n")
        rep.write("FIN DEL REPORTE\n")

    # ====================================================
    # FIGURA 3: MATRIZ DE ROBUSTEZ (GAPS)
    # ====================================================
    print("Generando Figura 3: Matriz de Robustez...")
    df_v10 = df[df['v_frag'] == 10]
    unique_rgap = sorted(df_v10['r_gap'].unique())
    unique_mgap = sorted(df_v10['M_gap'].unique())
    
    if len(unique_rgap) > 0 and len(unique_mgap) > 0:
        fig, axes = plt.subplots(len(unique_mgap), len(unique_rgap), 
                                 figsize=(4 * len(unique_rgap), 4 * len(unique_mgap)), 
                                 sharex=True, sharey=True)
        
        # Asegurar 2D
        axes = np.atleast_2d(axes)
        if len(unique_mgap) == 1: axes = axes.T if len(unique_rgap) > 1 else axes
        
        for i, mg in enumerate(unique_mgap):
            for j, rg in enumerate(unique_rgap):
                ax = axes[i, j]
                df_panel = df_v10[(df_v10['M_gap'] == mg) & (df_v10['r_gap'] == rg)]
                if len(df_panel) > 0:
                    # Usamos la misma función de binning pero quitamos el título global y ajustamos tamaño
                    sns.scatterplot(
                        data=df_panel, x='log10_Memb0', y='frac_h2o_percent',
                        hue='alpha', palette='viridis', alpha=0.3, ax=ax,
                        edgecolor=None, s=20, legend=False
                    )
                    medians = df_panel.groupby(['alpha', 'log10_Memb0'])['frac_h2o_percent'].median().reset_index()
                    sns.lineplot(
                        data=medians, x='log10_Memb0', y='frac_h2o_percent',
                        hue='alpha', palette='viridis', ax=ax, marker='o', 
                        linewidth=1.5, markersize=6, legend=False
                    )
                
                ax.set_title(f"$M_{{gap}}={mg}, R_{{gap}}={rg}$", fontsize=10)
                if i == len(unique_mgap) - 1:
                    ax.set_xlabel(r"$\log_{10}(M_{\rm emb, 0}/M_\oplus)$", fontsize=10)
                else:
                    ax.set_xlabel("")
                    
                if j == 0:
                    ax.set_ylabel(r"$f_{\rm H_2O}$ (%)", fontsize=10)
                else:
                    ax.set_ylabel("")
                    
                ax.tick_params(direction="in", top=True, right=True)

        fig.tight_layout()
        fig.savefig("data/benchmarks/scatter_matriz_gap_agua.png", dpi=300)
        plt.close(fig)
    else:
        print("Advertencia: No hay suficientes datos de gaps (M_gap, r_gap) para la matriz.")

    print("\n¡Pipeline estadístico completado exitosamente!")
    print("Revisar: 'data/benchmarks/' y 'data/estadisticas_agua_reporte.txt'")

if __name__ == '__main__':
    main()

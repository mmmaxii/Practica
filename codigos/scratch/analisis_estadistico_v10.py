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
    sns.scatterplot(
        data=df, x='log10_Memb0', y='frac_h2o_percent',
        hue='alpha', palette='viridis', alpha=0.3, ax=ax,
        edgecolor=None, s=40, legend=False
    )
    medians = df.groupby(['alpha', 'log10_Memb0'])['frac_h2o_percent'].median().reset_index()
    sns.lineplot(
        data=medians, x='log10_Memb0', y='frac_h2o_percent',
        hue='alpha', palette='viridis', ax=ax, marker='o', 
        linewidth=2.5, markersize=8, legend='full'
    )
    ax.set_title(title)
    ax.set_xlabel(r"$\log_{10}(M_{\rm emb, 0}/M_\oplus)$")
    ax.set_ylabel(r"Fracción de Agua Final (%)")
    handles, labels = ax.get_legend_handles_labels()
    if len(handles) > 0:
        ax.legend(handles, labels, title=r"$\alpha$", loc='upper right', bbox_to_anchor=(1.0, 1.0))

def main():
    print("Cargando datos y aislando v_frag = 10 m/s...")
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print("Error: No se encontró el CSV.")
        return
        
    df = pd.read_csv(csv_path)
    # Filtro Exclusivo
    df = df[(df['M_emb0'] > 0) & (df['alpha'] > 0) & (df['v_frag'] == 10)].copy()
    
    df['log10_Memb0'] = np.log10(df['M_emb0'])
    df['log10_alpha'] = np.log10(df['alpha'])
    
    report_path = "data/estadisticas_agua_reporte_v10.txt"
    with open(report_path, "w", encoding='utf-8') as rep:
        rep.write("==================================================\n")
        rep.write(" REPORTE AISLADO: SOLO v_frag = 10 m/s\n")
        rep.write("==================================================\n\n")
        rep.write(f"Total de simulaciones procesadas: {len(df)}\n\n")
        
        # 1. Correlaciones Simples
        p_corr, p_pval = pearsonr(df['log10_Memb0'], df['frac_h2o_percent'])
        s_corr, s_pval = spearmanr(df['log10_Memb0'], df['frac_h2o_percent'])
        
        rep.write("--- 1. CORRELACIONES BIVARIADAS ---\n")
        rep.write(f"Correlación Pearson (Lineal):  r = {p_corr:.4f} (p-value = {p_pval:.2e})\n")
        rep.write(f"Correlación Spearman (Monótona): rho = {s_corr:.4f} (p-value = {s_pval:.2e})\n\n")

        # Figura: Scatter Exclusivo
        print("Generando Figura Scatter Exclusiva...")
        fig, ax = plt.subplots(figsize=(8, 6))
        plot_scatter_binned(df, ax, "Correlación Pura (Aislado a $v_{\\rm frag} = 10$ m/s)")
        fig.tight_layout()
        fig.savefig("data/benchmarks/scatter_global_agua_v10.png", dpi=300)
        plt.close(fig)

        # Regresiones
        print("Calculando OLS...")
        scaler = StandardScaler()
        # Excluimos v_frag porque ahora es constante y generaría colinealidad perfecta (singular matrix)
        vars_to_scale = ['log10_Memb0', 'log10_alpha', 'M_gap', 'r_gap']
        df_scaled = df.copy()
        df_scaled[vars_to_scale] = scaler.fit_transform(df[vars_to_scale])
        
        model_simple = smf.ols("frac_h2o_percent ~ log10_Memb0 + log10_alpha + M_gap + r_gap", data=df_scaled).fit()
        model_ext = smf.ols("frac_h2o_percent ~ log10_Memb0 * log10_alpha + M_gap + r_gap", data=df_scaled).fit()
        
        rep.write("--- 2. REGRESIÓN MÚLTIPLE OLS (Variables Estandarizadas) ---\n")
        
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

        # Random Forest
        print("Entrenando Random Forest...")
        X = df_scaled[vars_to_scale]
        y = df['frac_h2o_percent']
        rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1)
        rf.fit(X, y)
        result = permutation_importance(rf, X, y, n_repeats=10, random_state=42, n_jobs=-1)
        
        rep.write("--- 3. IMPORTANCIA RELATIVA (Random Forest Permutation Importance) ---\n")
        imp_df = pd.DataFrame({'Variable': vars_to_scale, 'Importancia': result.importances_mean})
        imp_df = imp_df.sort_values(by='Importancia', ascending=False)
        for idx, row in imp_df.iterrows():
            rep.write(f"{row['Variable']:>15}: {row['Importancia']:.4f}\n")
            
        rep.write("\n==================================================\n")
        rep.write("FIN DEL REPORTE V10\n")

    print("\n¡Análisis Exclusivo v10 completado!")
    print("Revisar: 'data/benchmarks/scatter_global_agua_v10.png' y 'data/estadisticas_agua_reporte_v10.txt'")

if __name__ == '__main__':
    main()

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier, export_text
from sklearn.model_selection import cross_val_score
from sklearn.metrics import f1_score, make_scorer

plt.rcParams.update({
    "text.usetex": False,
    "font.family": "serif",
    "axes.labelsize": 12,
    "font.size": 12,
    "xtick.direction": "in",
    "ytick.direction": "in",
})

def main():
    print("Iniciando Machine Learning - Etapa 2: Heatmaps, Árboles y Ablación")
    csv_path = "data/summary_master_todos_casos.csv"
    df = pd.read_csv(csv_path)
    
    df = df[(df['M_emb0'] > 0) & (df['alpha'] > 0)].copy()
    df['log10_Memb0'] = np.log10(df['M_emb0'])
    df['log10_alpha'] = np.log10(df['alpha'])
    
    df['waterworld'] = ((df['M_final'] > 0.1) & 
                        (df['frac_h2o_percent'] > 10) & 
                        (df['frac_h2o_percent'] < 20)).astype(int)
                        
    report_path = "data/ml_stage2_report.txt"
    with open(report_path, "w", encoding='utf-8') as f:
        f.write("==================================================\n")
        f.write("   MACHINE LEARNING STAGE 2: DESCRUBRIMIENTO\n")
        f.write("==================================================\n\n")
        
        # 1. ¿Dónde están los 90 Waterworlds?
        ww_df = df[df['waterworld'] == 1]
        f.write("--- DISTRIBUCIÓN DE LOS WATERWORLDS EN EL DATASET ---\n")
        f.write("Distribución por alpha:\n")
        f.write(ww_df.groupby('alpha').size().to_string() + "\n\n")
        f.write("Distribución por M_gap:\n")
        f.write(ww_df.groupby('M_gap').size().to_string() + "\n\n")
        f.write("Distribución por log10_Memb0:\n")
        f.write(ww_df.groupby('log10_Memb0').size().to_string() + "\n\n")
        
        # 2. Heatmaps de Éxito
        print("Generando Heatmaps...")
        # Pivot alpha vs v_frag
        pivot_vfrag = df.pivot_table(values='waterworld', index='alpha', columns='v_frag', aggfunc='mean') * 100
        plt.figure(figsize=(6, 5))
        sns.heatmap(pivot_vfrag, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': 'P(Waterworld) %'})
        plt.title("Éxito de Waterworld: $\\alpha$ vs $v_{\\rm frag}$")
        plt.tight_layout()
        plt.savefig("data/benchmarks/heatmap_alpha_vfrag.png", dpi=300)
        plt.close()
        
        # Pivot alpha vs M_gap (TODOS los v_frag)
        pivot_mgap = df.pivot_table(values='waterworld', index='alpha', columns='M_gap', aggfunc='mean') * 100
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_mgap, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': 'P(Waterworld) %'})
        plt.title("Éxito de Waterworld: $\\alpha$ vs $M_{\\rm gap}$ (Todos los $v_{\\rm frag}$)")
        plt.tight_layout()
        plt.savefig("data/benchmarks/heatmap_alpha_mgap.png", dpi=300)
        plt.close()
        
        # Pivot alpha vs M_gap (SOLO v_frag = 10)
        df_v10 = df[df['v_frag'] == 10]
        pivot_mgap_v10 = df_v10.pivot_table(values='waterworld', index='alpha', columns='M_gap', aggfunc='mean') * 100
        plt.figure(figsize=(8, 6))
        sns.heatmap(pivot_mgap_v10, annot=True, fmt=".1f", cmap="YlGnBu", cbar_kws={'label': 'P(Waterworld) %'})
        plt.title("Éxito de Waterworld: $\\alpha$ vs $M_{\\rm gap}$ ($v_{\\rm frag} = 10$ m/s)")
        plt.tight_layout()
        plt.savefig("data/benchmarks/heatmap_alpha_mgap_v10.png", dpi=300)
        plt.close()

        # 3. Decision Tree Classifier
        print("Entrenando Decision Tree...")
        features = ['log10_alpha', 'v_frag', 'M_gap', 'r_gap', 'log10_Memb0']
        X = df[features]
        y = df['waterworld']
        
        dt = DecisionTreeClassifier(max_depth=3, class_weight='balanced', random_state=42)
        dt.fit(X, y)
        
        f.write("--- ÁRBOL DE DECISIÓN (Reglas de Formación) ---\n")
        f.write("Max Depth = 3, Class Weight = Balanced\n\n")
        rules = export_text(dt, feature_names=features)
        f.write(rules)
        f.write("\n")
        
        # 4. Ablation Study
        print("Realizando Ablation Study...")
        rf = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42, n_jobs=-1)
        
        # Función para evaluar F1
        def eval_features(feat_list):
            scores = cross_val_score(rf, df[feat_list], y, cv=5, scoring='f1')
            return scores.mean()
            
        f1_full = eval_features(features)
        
        feat_no_memb = [c for c in features if c != 'log10_Memb0']
        f1_no_memb = eval_features(feat_no_memb)
        
        feat_no_alpha = [c for c in features if c != 'log10_alpha']
        f1_no_alpha = eval_features(feat_no_alpha)
        
        feat_no_mgap = [c for c in features if c != 'M_gap']
        f1_no_mgap = eval_features(feat_no_mgap)
        
        f.write("--- ESTUDIO DE ABLACIÓN (Importancia Causal) ---\n")
        f.write(f"Métrica: F1-Score (Cross-Validation K=5)\n\n")
        f.write(f"Modelo COMPLETO (Todas las variables): F1 = {f1_full:.4f}\n")
        f.write(f"Modelo SIN Memb0 (Se quita embrión)  : F1 = {f1_no_memb:.4f} (Caída: {f1_full - f1_no_memb:.4f})\n")
        f.write(f"Modelo SIN alpha (Se quita turb.)    : F1 = {f1_no_alpha:.4f} (Caída: {f1_full - f1_no_alpha:.4f})\n")
        f.write(f"Modelo SIN M_gap (Se quita planeta)  : F1 = {f1_no_mgap:.4f} (Caída: {f1_full - f1_no_mgap:.4f})\n")
        
        f.write("\nInterpretación: Si la caída es masiva, el modelo colapsa sin esa variable. Si la caída es cero, la variable no aportaba información única.\n")

    print("\n¡Etapa 2 Completada! Archivos generados:")
    print("- data/benchmarks/heatmap_alpha_vfrag.png")
    print("- data/benchmarks/heatmap_alpha_mgap.png")
    print("- data/ml_stage2_report.txt")

if __name__ == '__main__':
    main()

import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier

def main():
    print("Iniciando Machine Learning - Etapa 1: Feature Importances Puras")
    csv_path = "data/summary_master_todos_casos.csv"
    if not os.path.exists(csv_path):
        print(f"Error: No se encontró {csv_path}")
        return
        
    df = pd.read_csv(csv_path)
    
    # Pre-procesamiento de features
    df = df[(df['M_emb0'] > 0) & (df['alpha'] > 0)].copy()
    df['log10_Memb0'] = np.log10(df['M_emb0'])
    df['log10_alpha'] = np.log10(df['alpha'])
    
    # Definición del Target C (Waterworld)
    df['waterworld'] = (df['M_final'] > 0.1) & (df['frac_h2o_percent'] > 10) & (df['frac_h2o_percent'] < 20)
    df['waterworld'] = df['waterworld'].astype(int)
    
    # Features
    features = ['log10_alpha', 'v_frag', 'M_gap', 'r_gap', 'log10_Memb0']
    X = df[features]
    
    # Target A: M_final
    y_A = df['M_final']
    rf_A = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    rf_A.fit(X, y_A)
    
    # Target B: frac_h2o_percent
    y_B = df['frac_h2o_percent']
    rf_B = RandomForestRegressor(n_estimators=200, random_state=42, n_jobs=-1)
    rf_B.fit(X, y_B)
    
    # Target C: waterworld
    y_C = df['waterworld']
    rf_C = RandomForestClassifier(n_estimators=200, random_state=42, n_jobs=-1)
    rf_C.fit(X, y_C)
    
    # Reportar
    report_path = "data/ml_stage1_importances.txt"
    with open(report_path, "w", encoding='utf-8') as f:
        f.write("==================================================\n")
        f.write("   MACHINE LEARNING STAGE 1: RANDOM FOREST\n")
        f.write("==================================================\n\n")
        f.write(f"Total simulaciones: {len(df)}\n")
        f.write(f"Total Waterworlds encontrados: {df['waterworld'].sum()} ({(df['waterworld'].sum() / len(df)) * 100:.1f}%)\n\n")
        
        def print_importances(model, target_name):
            f.write(f"--- OBJETIVO: {target_name} ---\n")
            importances = model.feature_importances_
            imp_df = pd.DataFrame({'Variable': features, 'Importancia': importances})
            imp_df = imp_df.sort_values(by='Importancia', ascending=False)
            for _, row in imp_df.iterrows():
                f.write(f"{row['Variable']:>15}: {row['Importancia']:.4f}\n")
            f.write("\n")
            
        print_importances(rf_A, "A. Masa Final (M_final) [Regressor]")
        print_importances(rf_B, "B. Fracción de Agua (f_H2O) [Regressor]")
        print_importances(rf_C, "C. Waterworld (Clasificador) [Classifier]")
        
    print(f"Completado. Reporte guardado en {report_path}.")

if __name__ == '__main__':
    main()

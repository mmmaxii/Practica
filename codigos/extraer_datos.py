import os
import csv
import glob
from pipeline_analysis.analyzer_base import GeneralAnalyzer

def compilar_csv_burbujas():
    base_dir = "data/runs"
    m0_earth = 0.01 # Masa inicial semilla ultra-baja (sin sobreestimación)
    
    dir_general = os.path.join(base_dir, "vf_10ms", "general")
    cache_files = glob.glob(os.path.join(dir_general, "cache_general_a_*_M0_0.01.pkl"))
    
    todos_los_datos = []
    
    for cache_path in cache_files:
        basename = os.path.basename(cache_path)
        parts = basename.split('_')
        alpha_target = float(parts[3])
        
        analyzer = GeneralAnalyzer(dir_general, cache_path, m0_earth, alpha_target)
        datos = analyzer.load_or_extract()
        for d in datos:
            d['alpha'] = alpha_target
        todos_los_datos.extend(datos)
        
    # Añadir los datos del caso 'rounded' (ahora dinámico para todos los alphas)
    dir_rounded = os.path.join(base_dir, "vf_10ms", "rounded")
    cache_rounded_files = glob.glob(os.path.join(dir_rounded, "cache_rounded_a_*_M0_0.01.pkl"))
    
    for cache_path in cache_rounded_files:
        basename = os.path.basename(cache_path)
        parts = basename.split('_')
        # El nombre es cache_rounded_a_0.001_M0_0.01.pkl -> partes: [cache, rounded, a, 0.001, M0, 0.01.pkl]
        alpha_target = float(parts[3])
        
        analyzer_round = GeneralAnalyzer(dir_rounded, cache_path, m0_earth, alpha_target)
        datos_round = analyzer_round.load_or_extract()
        for d in datos_round:
            d['alpha'] = alpha_target
        todos_los_datos.extend(datos_round)
        
    # Escribir el CSV consolidado
    csv_out = os.path.join(base_dir, "summary_all_alphas_v10.csv")
    
    print(f"Escribiendo {len(todos_los_datos)} filas en {csv_out}...")
    with open(csv_out, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['run_name', 'alpha', 'r_gap', 'M_gap', 'M_final', 'frac_h2o'])
        for row in todos_los_datos:
            run_name = row['run_name']
            r_gap = row['r_gap']
            m_gap = row['M_gap']
            alpha = row.get('alpha', 0.001)
            m_final = row['mass_e'][-1] if len(row['mass_e']) > 0 else 0.0
            frac_h2o = row['frac_h2o_final'] / 100.0 
            
            writer.writerow([run_name, alpha, r_gap, m_gap, m_final, frac_h2o])
            
    print("Extracción y consolidación CSV terminada.")

if __name__ == "__main__":
    compilar_csv_burbujas()

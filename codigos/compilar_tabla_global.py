import os
import glob
import pickle
import csv
import re

def compilar_tabla_maestra():
    base_dir = "data/runs"
    # Buscar todos los caches dentro de las carpetas vf_*ms
    cache_files = glob.glob(os.path.join(base_dir, "vf_*ms", "**", "cache_*.pkl"), recursive=True)
    
    out_csv = os.path.join("data", "summary_master_todos_casos.csv")
    todas_las_filas = []
    
    print(f"Encontrados {len(cache_files)} archivos de caché para procesar.")
    
    for cache_path in cache_files:
        # Extraer información del path y del nombre del archivo
        # path típico: data/runs/vf_10ms/general/cache_general_a_0.001_M0_0.01.pkl
        path_parts = cache_path.replace("\\", "/").split('/')
        
        # Extraer v_frag de la carpeta 'vf_10ms'
        v_frag_folder = [p for p in path_parts if p.startswith('vf_')][0]
        v_frag_match = re.search(r'vf_(\d+)ms', v_frag_folder)
        v_frag = int(v_frag_match.group(1)) if v_frag_match else 10
        
        # Extraer escenario (general, rounded, etc)
        escenario = path_parts[path_parts.index(v_frag_folder) + 1]
        
        basename = os.path.basename(cache_path)
        # Regex para sacar alpha y M0 del nombre
        # ej: cache_general_a_0.001_M0_0.01.pkl
        match = re.search(r'_a_([0-9.]+)_M0_([0-9.eE-]+)\.pkl', basename)
        if match:
            alpha = float(match.group(1))
            m0 = float(match.group(2))
        else:
            print(f"No se pudo parsear {basename}. Saltando...")
            continue
            
        with open(cache_path, 'rb') as f:
            data = pickle.load(f)
            
        for d in data:
            run_name = d.get('run_name', 'unknown')
            r_gap = d.get('r_gap', 0.0)
            m_gap = d.get('M_gap', 0.0)
            
            # Masa final (último valor del array mass_e)
            mass_arr = d.get('mass_e', [])
            m_final = mass_arr[-1] if len(mass_arr) > 0 else m0
            
            # Parametros sinusoidales
            n_gaps = d.get('n_gaps', 0)
            amp_val = d.get('amp_val', 0.0)
            
            # Fracción de agua
            frac_h2o = d.get('frac_h2o_final', 0.0)
            
            todas_las_filas.append({
                'run_name': run_name,
                'v_frag': v_frag,
                'scenario': escenario,
                'alpha': alpha,
                'M_emb0': m0,
                'r_gap': r_gap,
                'M_gap': m_gap,
                'n_gaps': n_gaps,
                'amp_val': amp_val,
                'M_final': m_final,
                'frac_h2o_percent': frac_h2o
            })
            
    print(f"Escribiendo {len(todas_las_filas)} entradas en {out_csv}...")
    with open(out_csv, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'run_name', 'v_frag', 'scenario', 'alpha', 'M_emb0', 
            'r_gap', 'M_gap', 'n_gaps', 'amp_val', 'M_final', 'frac_h2o_percent'
        ])
        writer.writeheader()
        writer.writerows(todas_las_filas)
        
    print("¡Tabla Maestra generada con éxito!")

if __name__ == '__main__':
    compilar_tabla_maestra()

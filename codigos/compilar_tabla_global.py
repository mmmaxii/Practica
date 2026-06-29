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
            
            mass_arr = d.get('mass_e', [])
            if len(mass_arr) < 2:
                m_final = m0
                m_prev = m0
            else:
                m_final = mass_arr[-1]
                m_prev = mass_arr[-2]
                
            n_gaps = d.get('n_gaps', 0)
            amp_val = d.get('amp_val', 0.0)
            frac_h2o = d.get('frac_h2o_final', 0.0)
            f_water_final = frac_h2o / 100.0
            
            t_final = d.get('times_yr', [0])[-1] / 1e6 if 'times_yr' in d else 0.0
            
            # Extracción robusta de M_iso_e
            m_iso = d.get('m_iso_e', 3.955)
            if m_iso is None:
                m_iso = 3.955
            elif isinstance(m_iso, (list, tuple)):
                m_iso = m_iso[-1] if len(m_iso) > 0 else 3.955
            elif type(m_iso).__name__ == 'ndarray':
                m_iso = m_iso[-1] if len(m_iso) > 0 else 3.955
                
            if m_iso >= 100.0:
                m_iso = 3.955
                
            R_var = abs(m_final - m_prev) / m_final if m_final > 0 else 0.0
            
            # Lógica de Validación
            is_valid = True
            discard_reason = ""
            
            if m_final <= 0:
                is_valid = False
                discard_reason = "Fell into star (m <= 0)"
                m_final = 0.0
            else:
                if m_final >= m_iso * 0.90:
                    pass # Si alcanza el 90% de la masa de aislamiento, es válido sin importar el tiempo
                elif t_final >= 5.0:
                    pass # Si duró 5 Myr, es válido
                elif t_final >= 1.0 and t_final < 5.0:
                    if R_var >= 0.10:
                        is_valid = False
                        discard_reason = "Unstable (R >= 10%)"
                else:
                    is_valid = False
                    discard_reason = "Premature (< 1 Myr and < 90% M_iso)"
                            
            # Categorización
            cat_masa = ""
            cat_agua = ""
            if m_final < 0.1:
                cat_masa = "< 0.1 M_E"
            elif m_final >= m_iso * 0.90:
                cat_masa = "Isolation Mass"
            else:
                cat_masa = "> 0.1 M_E (No iso)"
                
            if cat_masa in ["Isolation Mass", "> 0.1 M_E (No iso)"]:
                if f_water_final == 0.0:
                    cat_agua = "0% Agua"
                elif 0.0 < f_water_final < 0.0001:
                    cat_agua = "Mundos Áridos"
                elif 0.0001 <= f_water_final < 0.001:
                    if m_final <= 1.0:
                        cat_agua = "Análogo Terrestre Estricto"
                    else:
                        cat_agua = "Súper-Tierra Seca"
                elif 0.001 <= f_water_final <= 0.10:
                    cat_agua = "Transición (0.1%-10%)"
                else:
                    cat_agua = "Water Worlds (>10%)"
                    
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
                'M_prev': m_prev,
                'T_final_Myr': t_final,
                'M_iso_E': m_iso,
                'R_variation': R_var,
                'frac_h2o_percent': frac_h2o,
                'Is_Valid': is_valid,
                'Discard_Reason': discard_reason,
                'Cat_Masa': cat_masa,
                'Cat_Agua': cat_agua
            })
            
    print(f"Escribiendo {len(todas_las_filas)} entradas en {out_csv}...")
    try:
        with open(out_csv, 'w', newline='', encoding='utf-8') as f:
            fieldnames = [
                'run_name', 'v_frag', 'scenario', 'alpha', 'M_emb0', 
                'r_gap', 'M_gap', 'n_gaps', 'amp_val', 'M_final', 'M_prev',
                'T_final_Myr', 'M_iso_E', 'R_variation', 'frac_h2o_percent',
                'Is_Valid', 'Discard_Reason', 'Cat_Masa', 'Cat_Agua'
            ]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(todas_las_filas)
        print("¡Tabla Maestra generada con éxito!")
    except Exception as e:
        print(f"Error escribiendo {out_csv}: {e}")

if __name__ == '__main__':
    compilar_tabla_maestra()

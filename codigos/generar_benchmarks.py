import os
import csv
import numpy as np
from pipeline_analysis.snapshot_analyzer import SnapshotAnalyzer
from pipeline_analysis.plotters import PlotterBenchmarks

# CONFIGURACIÓN DE EJECUCIÓN
MODO_MOCK = False  
DIRECTORIO_SALIDA = "data/benchmarks/"
os.makedirs(DIRECTORIO_SALIDA, exist_ok=True)

print("=== INICIANDO PIPELINE DE BENCHMARKS DEFINITIVOS ===")

# ==========================================
# STEP 1: GENERACIÓN DEL BUBBLE CHART (ALPHA = 0.001)
# ==========================================
print("\n[1/3] Generando Bubble Chart para alpha = 0.001...")

csv_path = "data/runs/summary_all_alphas_v10.csv"
data_por_alpha = {}

if os.path.exists(csv_path):
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            alpha = row['alpha']
            if alpha not in data_por_alpha:
                data_por_alpha[alpha] = {'r_gap': [], 'M_gap': [], 'M_final': [], 'frac_h2o': []}
                
            data_por_alpha[alpha]['r_gap'].append(float(row['r_gap']))
            data_por_alpha[alpha]['M_gap'].append(float(row['M_gap']))
            data_por_alpha[alpha]['M_final'].append(float(row['M_final']))
            data_por_alpha[alpha]['frac_h2o'].append(float(row['frac_h2o']))
            
    df_por_alpha = {}
    for alpha_str, arrays in data_por_alpha.items():
        df_bubble = {
            'r_gap': np.array(arrays['r_gap']),
            'M_gap': np.array(arrays['M_gap']),
            'M_final': np.array(arrays['M_final']),
            'frac_h2o': np.array(arrays['frac_h2o'])
        }
        
        # Validar si df_bubble tiene datos
        if len(df_bubble['M_final']) == 0: continue
        
        df_por_alpha[alpha_str] = df_bubble
        
        PlotterBenchmarks.plot_bubble_chart(
            df_bubble, 
            os.path.join(DIRECTORIO_SALIDA, f"bubble_chart_alpha_{alpha_str}.png"), 
            alpha_val=f"{alpha_str}"
        )
        print(f"-> Bubble Chart (alpha={alpha_str}) guardado con éxito.")
        
    # Generar el Mosaico estilo Paper (para todos los alphas)
    # Ordenar los alphas numéricamente para que el mosaico tenga sentido físico
    alphas_mosaic = sorted(list(df_por_alpha.keys()), key=float)
    if alphas_mosaic:
        PlotterBenchmarks.plot_bubble_chart_mosaic(
            df_por_alpha,
            alphas_mosaic,
            os.path.join(DIRECTORIO_SALIDA, "bubble_chart_mosaico_paper.png")
        )
        print(f"-> Mosaico de Paper generado exitosamente para todos los alphas: {alphas_mosaic}.")
else:
    print(f"-> ADVERTENCIA: No se encontró {csv_path}. Ejecuta extraer_datos.py primero.")

# ==========================================
# STEP 2 & 3: DIAGRAMAS DE HOVMÖLLER Y PEBBLE FLUX
# ==========================================
print("\n[2/3] Procesando hidrodinámica de Casos Arquetípicos...")

arquetipos = [
    {
        "path": "data/runs/vf_3ms/general/run_r15.0_m2.0_a0.001_v3.0",
        "name": "Sweet Spot",
        "title": "Sweet Spot (M_gap=2.0, R_gap=15, v_frag=3)"
    },
    {
        "path": "data/runs/vf_1ms/general/run_r15.0_m2.0_a0.001_v1.0",
        "name": "Fracaso",
        "title": "Asfixia (M_gap=2.0, R_gap=15, v_frag=1)"
    },
    {
        "path": "data/runs/vf_10ms/general/run_r15.0_m0.1_a0.001",
        "name": "Exceso",
        "title": "Gigante Helado (M_gap=0.1, R_gap=15, v_frag=10)"
    }
]

for arq in arquetipos:
    print(f"-> Analizando Arquetipo: {arq['name']}...")
    if os.path.exists(arq['path']):
        analyzer = SnapshotAnalyzer(arq['path'])
        t_arr, r_arr, eps_mat, flux_arr, amax_arr = analyzer.extract_spatiotemporal_data(subsampling=1)
        
        # Aproximación analítica de Oka para R_ice(t)
        r_ice_real = 2.73 * (t_arr / (0.2 * 3.15576e7 * 1e6))**(-0.5) 
        r_ice_real = np.clip(r_ice_real, 0.5, 2.73)
        
        PlotterBenchmarks.plot_hovmoller(
            t_arr, r_arr, eps_mat, r_ice_real, 
            os.path.join(DIRECTORIO_SALIDA, f"hovmoller_{arq['name'].replace(' ', '_').lower()}.png"), 
            arq['title']
        )
        PlotterBenchmarks.plot_pebble_flux(
            t_arr, flux_arr, 
            os.path.join(DIRECTORIO_SALIDA, f"pebble_flux_{arq['name'].replace(' ', '_').lower()}.png"), 
            arq['title']
        )
    else:
        print(f"   [!] Directorio no encontrado: {arq['path']}")

# ==========================================
# STEP 4: MOSAICOS 3x2 DE HOVMÖLLER Y PEBBLE FLUX
# ==========================================
print("\n[3/3] Construyendo Mosaicos Spatio-Temporales (Hovmöller y Flux) para todos los alphas...")

data_hov = {}
for alpha in alphas_mosaic:
    arq_path = f"data/runs/vf_10ms/general/run_r10.0_m0.3_a{alpha}"
    if os.path.exists(arq_path):
        analyzer = SnapshotAnalyzer(arq_path)
        t_arr, r_arr, eps_mat, flux_arr, amax_arr = analyzer.extract_spatiotemporal_data(subsampling=1)
        data_hov[alpha] = {
            't_arr': t_arr,
            'r_arr': r_arr,
            'eps_mat': eps_mat,
            'flux_arr': flux_arr,
            'amax_arr': amax_arr
        }
    else:
        print(f"   [!] Faltan datos espaciotemporales para alpha={alpha}: {arq_path}")

if data_hov:
    PlotterBenchmarks.plot_hovmoller_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "hovmoller_mosaico_alphas.png")
    )
    print("-> Mosaico de Hovmöller generado exitosamente.")
    
    PlotterBenchmarks.plot_pebble_flux_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "pebble_flux_mosaico_alphas.png")
    )
    print("-> Mosaico de Pebble Flux (solo flujo) generado exitosamente.")
    
    PlotterBenchmarks.plot_flux_amax_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "pebble_flux_amax_mosaico.png")
    )
    print("-> Mosaico de Flujo + A_max generado exitosamente.")

print("\n=== PIPELINE EJECUTADO: BENCHMARKS GUARDADOS EN:", DIRECTORIO_SALIDA, "===")

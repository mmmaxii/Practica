import os
from pipeline_analysis.snapshot_analyzer import SnapshotAnalyzer
from pipeline_analysis.plotters import PlotterBenchmarks

DIRECTORIO_SALIDA = "data/benchmarks/"
os.makedirs(DIRECTORIO_SALIDA, exist_ok=True)

print("=== INICIANDO PIPELINE DE MOSAICOS PARA v_frag = 3 m/s ===")

alphas_mosaic = ['0.0001', '0.0003', '0.0005', '0.0007', '0.001']
data_hov = {}

print("Extrayendo datos espaciotemporales (Hovmöller y Flux) para v_frag = 3 m/s...")

for alpha in alphas_mosaic:
    arq_path = f"data/runs/vf_3ms/general/run_r10.0_m0.3_a{alpha}_v3.0"
    if os.path.exists(arq_path):
        print(f" -> Procesando: {arq_path}")
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
    print("\nGenerando mosaico de Hovmöller...")
    PlotterBenchmarks.plot_hovmoller_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "hovmoller_mosaico_alphas_v3.png")
    )
    
    print("Generando mosaico de Pebble Flux...")
    PlotterBenchmarks.plot_pebble_flux_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "pebble_flux_mosaico_alphas_v3.png")
    )
    
    print("Generando mosaico de Pebble Flux + A_max...")
    PlotterBenchmarks.plot_flux_amax_mosaic(
        data_hov, 
        alphas_mosaic, 
        os.path.join(DIRECTORIO_SALIDA, "pebble_flux_amax_mosaico_v3.png")
    )
    
    print("\n[OK] Mosaicos para v_frag = 3 m/s generados exitosamente en data/benchmarks/")
else:
    print("\n[ERROR] No se pudo encontrar ningún dato para generar los mosaicos.")

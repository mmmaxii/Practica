# -*- coding: utf-8 -*-
import os
import sys
import numpy as np

# Asegurar que importamos PA3Py desde la raíz
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
import dustpy.constants as c

def analizar_casos(base_dir="runs_geryon/benchmarks"):
    # Rutas esperadas de las carpetas generadas
    caso_05 = os.path.join(base_dir, "run_ngap7_A1.0_a0.001_rmin0.5")
    caso_07 = os.path.join(base_dir, "run_ngap7_A1.0_a0.001_rmin0.7")
    
    m0_earth = 0.01  # Masa inicial de semilla típica
    
    print("==================================================")
    print("Analizando Benchmark A/B Test de r_min")
    print("==================================================")
    
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    
    for nombre, ruta, color in [("Caso A (0.5 AU)", caso_05, 'blue'), ("Caso B (0.7 AU)", caso_07, 'orange')]:
        if not os.path.exists(os.path.join(ruta, "data0000.hdf5")):
            print(f"{nombre}: NO ENCONTRADO en {ruta} (O no ha terminado de descargar)")
            continue
            
        print(f"\nExtrayendo datos de {nombre}...")
        try:
            import contextlib
            with open(os.devnull, 'w') as fnull, contextlib.redirect_stdout(fnull):
                pa3 = PebbleAccretionModule3.from_datadir(ruta, M_star=1.0)
                res = pa3.run_growth([1.0], M0_g=m0_earth * c.M_earth)
                
            hist = res[1.0]
            if len(hist) > 0:
                t_yr = hist[:, 0] / c.year
                mass_e = hist[:, 1] / c.M_earth
                masa_final = mass_e[-1]
                print(f"-> Masa final del embrión a 1 AU: {masa_final:.4f} M_Earth")
                
                plt.plot(t_yr, mass_e, label=f"{nombre} (Final: {masa_final:.2f} M_E)", color=color, linewidth=2)
            else:
                print("-> Error: Historial vacío. El embrión no creció o la data está corrupta.")
        except Exception as e:
            print(f"-> Error procesando {nombre}: {e}")

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Tiempo [Años]')
    plt.ylabel(r'Masa del Embrión a 1 AU [$M_\oplus$]')
    plt.title('Comparación Benchmark r_min (N=7, A=1.0, a=0.001, v_frag=10 m/s)')
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend()
    
    img_path = os.path.join(base_dir, "benchmark_rmin_comparison.png")
    plt.savefig(img_path, dpi=300, bbox_inches='tight')
    print(f"\n-> Gráfico de evolución temporal guardado en: {img_path}")

if __name__ == "__main__":
    # Permite pasar la ruta de benchmarks por consola, si no usa la por defecto
    directorio = sys.argv[1] if len(sys.argv) > 1 else "runs_geryon/benchmarks"
    analizar_casos(directorio)

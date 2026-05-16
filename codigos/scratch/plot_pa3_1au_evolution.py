import os
import sys
import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors

# Configurar importación de PA3Py
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(_HERE), "PA3Py"))
from PebbleAccretion3 import PebbleAccretionModule3
import dustpy.constants as c

DATAROOT = "data/runs/1myr"
SAVEDIR = "data/figures/1myr/embryo_1au"
os.makedirs(SAVEDIR, exist_ok=True)

# Configuraciones dadas por el usuario
EMBRYOS_AU = [1]
CATEGORIES = ["sup_earth", "nep", "sat", "jup", "sup_jup"]
# Gaps mayores o iguales a 3.0
GAPS_AU = [3.0, 5.0, 7.0, 10.0, 15.0, 20.0]
M_INIT_EARTH = 0.1
T_MIN_YR = 100.0

def main():
    # Diccionario para almacenar la masa de núcleo (M_core) en cada t_yr
    # Esto servirá para calcular dM y guardarlo en Pandas
    # Esto servirá para calcular dM y guardarlo en Pandas
    history_dict = {}
    miso_dict = {} # <-- Diccionario para guardar M_iso
    time_dict = {} # <-- Guardamos t_yr por simulacion para evitar mismatch

    # Recolectar datos
    print("=======================================================")
    print(" 1. Extrayendo datos de simulaciones (Embrion @ 2.5 AU)")
    print("=======================================================")
    
    for cat in CATEGORIES:
        for gap in GAPS_AU:
            run_name = f"single_{cat}_{gap}au"
            run_path = os.path.join(DATAROOT, run_name)
            
            # Verificar que exista y tenga snapshots HDF5
            if not os.path.isdir(run_path) or not glob.glob(os.path.join(run_path, "data*.hdf5")):
                continue
                
            print(f"  → Procesando: {run_name}")
            try:
                pa3 = PebbleAccretionModule3.from_datadir(run_path, t_min_yr=T_MIN_YR)
            except Exception as e:
                print(f"    [!] Error al cargar: {e}")
                continue
            
            # Correr el crecimiento para 2.5 AU con la masa inicial indicada (pasada en gramos)
            results = pa3.run_growth(EMBRYOS_AU, M0_g=M_INIT_EARTH * c.M_earth)
            
            hist = results.get(1.0)
            if hist is None or len(hist) == 0:
                print(f"    [!] No hay historia para 1 AU en {run_name}")
                continue
            
            t_yr = hist[:, 0] / c.year  # Convertir de segundos a años
            M_core = hist[:, 1] / c.M_earth  # en masas terrestres
            
            # --- COMO LO HICE ---
            # El array hist devuelto por run_growth() tiene las siguientes columnas:
            # col 0: tiempo, col 1: M_core, col 2-4: masas volatiles/silicatos
            # col 5: r_snow, col 6: M_iso.
            # Por ende, extraemos M_iso tomando la columna 6 y dividiendo por M_earth.
            M_iso_arr = hist[:, 6] / c.M_earth
            
            # Guardamos la serie de masa acumulada por ahora
            history_dict[run_name] = M_core
            miso_dict[run_name] = M_iso_arr
            time_dict[run_name] = t_yr
                
    if not history_dict:
        print("No se encontraron simulaciones válidas.")
        return

    # El usuario pide guardar la MASA ACRETADA EN ESE INTERVALO (dM)
    # dM_i = M_i - M_{i-1}. Para i=0, dM_0 = M_0 - M_init
    
    # Lista de nombres de run encontrados en el orden que los vamos a escribir
    run_names = list(history_dict.keys())
    
    # Encontrar la simulación con más snapshots para usar como tiempo maestro del CSV
    master_time = max(time_dict.values(), key=len)
    
    csv_path = os.path.join(SAVEDIR, "accreted_mass_intervals_2.5au.csv")
    with open(csv_path, 'w') as f:
        # Escribir cabeceras
        header = "Time[yr]," + ",".join(run_names) + "\n"
        f.write(header)
        
        # Iterar sobre cada fila (snapshot)
        for i in range(len(master_time)):
            row = [str(master_time[i])]
            for run_name in run_names:
                M_core_arr = history_dict[run_name]
                if i < len(M_core_arr):
                    M_current = M_core_arr[i]
                    M_prev = M_core_arr[i-1] if i > 0 else M_INIT_EARTH
                    dM = M_current - M_prev
                    row.append(str(dM))
                else:
                    row.append("") # En caso de que haya una simulación con menos snapshots
            f.write(",".join(row) + "\n")
            
    print(f"\n[OK] Datos de intervalos (dM) guardados en: {csv_path}")

    # =================================================================
    # 2. Generar Gráficos Tipo 1: 1 plot por planeta, curvas por gap
    # =================================================================
    print("\nGenerando gráficos Tipo 1 (por planeta)...")
    cmap_gap = cm.get_cmap("viridis")
    norm_gap = mcolors.Normalize(vmin=min(GAPS_AU), vmax=max(GAPS_AU))
    
    for cat in CATEGORIES:
        plt.figure(figsize=(8, 6))
        plot_added = False
        
        for gap in GAPS_AU:
            run_name = f"single_{cat}_{gap}au"
            if run_name in history_dict:
                y_vals = history_dict[run_name]
                iso_vals = miso_dict[run_name]
                x_vals = time_dict[run_name]
                plt.plot(x_vals, y_vals, marker='o', markersize=3, 
                         color=cmap_gap(norm_gap(gap)), lw=2, label=f"Gap {gap} AU")
                # Graficar la masa de aislamiento con una linea punteada
                plt.plot(x_vals, iso_vals, linestyle='--', color=cmap_gap(norm_gap(gap)), lw=1.5, alpha=0.7)
                plot_added = True
                
        if plot_added:
            plt.xscale('log')
            plt.xlabel('Tiempo [años]', fontsize=12)
            plt.ylabel(r'Masa del núcleo $M_{core}$ [$M_{\oplus}$]', fontsize=12)
            plt.title(f'Evolución del embrión a 2.5 AU - Planeta: {cat.upper()}', fontsize=14)
            plt.grid(True, which="both", ls="--", alpha=0.4)
            plt.legend()
            
            fig_path = os.path.join(SAVEDIR, f"evol_2.5au_planeta_{cat}.png")
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()
            
    # =================================================================
    # 3. Generar Gráficos Tipo 2: 1 plot por gap, curvas por planeta
    # =================================================================
    print("Generando gráficos Tipo 2 (por distancia de gap)...")
    cmap_cat = cm.get_cmap("plasma")
    norm_cat = mcolors.Normalize(vmin=0, vmax=len(CATEGORIES)-1)
    
    for gap in GAPS_AU:
        plt.figure(figsize=(8, 6))
        plot_added = False
        
        for i, cat in enumerate(CATEGORIES):
            run_name = f"single_{cat}_{gap}au"
            if run_name in history_dict:
                y_vals = history_dict[run_name]
                iso_vals = miso_dict[run_name]
                x_vals = time_dict[run_name]
                plt.plot(x_vals, y_vals, marker='s', markersize=3, 
                         color=cmap_cat(norm_cat(i)), lw=2, label=f"Planeta: {cat}")
                # Graficar la masa de aislamiento con una linea punteada
                plt.plot(x_vals, iso_vals, linestyle='--', color=cmap_cat(norm_cat(i)), lw=1.5, alpha=0.7)
                plot_added = True
                
        if plot_added:
            plt.xscale('log')
            plt.xlabel('Tiempo [años]', fontsize=12)
            plt.ylabel(r'Masa del núcleo $M_{core}$ [$M_{\oplus}$]', fontsize=12)
            plt.title(f'Evolución del embrión a 2.5 AU - Gap a {gap} AU', fontsize=14)
            plt.grid(True, which="both", ls="--", alpha=0.4)
            plt.legend()
            
            fig_path = os.path.join(SAVEDIR, f"evol_2.5au_gap_{gap}au.png")
            plt.savefig(fig_path, dpi=300, bbox_inches='tight')
            plt.close()

    print("\n[OK] Todos los gráficos generados exitosamente.")

if __name__ == "__main__":
    main()

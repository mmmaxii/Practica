import numpy as np
import matplotlib.pyplot as plt
import dustpy.constants as c
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

def main():
    # Directorio de los datos
    datadir = "data/runs/5myr/single_nep_5.0au"
    
    print(f"Cargando datos desde {datadir}...")
    try:
        pa3 = PebbleAccretionModule3.from_datadir(datadir, M_star=1.0)
    except FileNotFoundError:
        print(f"Error: No se encontró el directorio {datadir}")
        return

    # Definir los 100 embriones desde 1.0 hasta 5.0 AU
    embryos_au = np.linspace(1.0, 5.0, 100)
    
    # Masa inicial: típicamente usamos 0.01 o 0.1 M_earth. Usaré 0.01 M_earth.
    m_init = 0.01 * c.M_earth
    
    print(f"Calculando crecimiento para {len(embryos_au)} embriones...")
    results = pa3.run_growth(embryos_au, M0_g=m_init)
    
    # Graficar
    plt.figure(figsize=(10, 6))
    
    distances = []
    final_masses = []
    iso_masses = []
    
    for r_au in sorted(results.keys()):
        hist = results[r_au]
        if len(hist) == 0:
            continue
            
        # Extraer masa final y masa de aislamiento del último snapshot
        # Columnas de hist: 0:t, 1:M_core, 2:H2O, 3:CO2, 4:silicates, 5:r_snow, 6:M_iso
        M_final = hist[-1, 1] / c.M_earth
        M_iso = hist[-1, 6] / c.M_earth
        
        distances.append(r_au)
        final_masses.append(M_final)
        iso_masses.append(M_iso)

    # Plot masa final
    plt.plot(distances, final_masses, 'o-', color='tab:blue', linewidth=2, label='Masa Final del Embrión')
    
    # Plot masa de aislamiento
    plt.plot(distances, iso_masses, '--', color='tab:red', linewidth=2, label='Masa de Aislamiento (M_iso)')

    plt.yscale('log')
    plt.xlabel('Distancia (AU)', fontsize=12)
    plt.ylabel(r'Masa ($M_\oplus$)', fontsize=12)
    plt.title('Masa Final de los Embriones vs Distancia (A 5 Myr)\nNeptuno a 5.0 AU', fontsize=14)
    
    # Limitar el eje Y al máximo de las masas finales de los embriones (como lo solicitó el usuario)
    max_mass = max(final_masses) if final_masses else 1.0
    plt.ylim(bottom=0.01, top=max_mass * 1.5)
    
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.legend(fontsize=12)
    
    # Guardar en la raíz
    output_file = "embryos_final_mass_vs_distance_nep.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Gráfico guardado exitosamente como '{output_file}' en la raíz del proyecto.")

if __name__ == "__main__":
    main()

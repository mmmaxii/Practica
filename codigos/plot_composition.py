import numpy as np
import matplotlib.pyplot as plt
import dustpy.constants as c
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

def main():
    # Usaremos los datos del Neptuno en 5.0 AU que acabamos de correr
    datadir = "data/runs/1myr/single_nep_10.0au"
    
    print(f"Cargando datos desde {datadir}...")
    try:
        pa3 = PebbleAccretionModule3.from_datadir(datadir, M_star=1.0)
    except FileNotFoundError:
        print(f"Error: No se encontró el directorio {datadir}")
        return

    # Estudiaremos un embrión que nazca en 2.0 AU
    r_emb = 2.0
    m_init = 0.01 * c.M_earth
    
    print(f"Calculando crecimiento y composición para embrión en {r_emb} AU...")
    results = pa3.run_growth([r_emb], M0_g=m_init)
    
    hist = results[r_emb]
    if len(hist) == 0:
        print("El embrión no registró historial.")
        return
        
    # Extraer columnas (ver docstring o PA3Py)
    # col 0: tiempo, col 1: M_core, col 2: H2O, col 3: CO2, col 4: silicates
    t_yr = hist[:, 0] / c.year
    M_total = hist[:, 1] / c.M_earth
    M_h2o = hist[:, 2] / c.M_earth
    M_silicates = hist[:, 4] / c.M_earth
    
    # Calcular fracción de hielo de agua (%)
    # Evitamos dividir por cero
    frac_h2o = np.where(M_total > 0, (M_h2o / M_total) * 100, 0.0)

    # ── Gráfico 1: Evolución de las Masas ──
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)
    
    ax1.plot(t_yr, M_total, lw=3, color='black', label='Masa Total')
    ax1.plot(t_yr, M_silicates, lw=2, color='saddlebrown', linestyle='--', label='Masa Silicatos (Roca)')
    ax1.plot(t_yr, M_h2o, lw=2, color='dodgerblue', linestyle='-.', label='Masa H2O (Hielo)')
    
    ax1.set_yscale('log')
    ax1.set_ylabel(r'Masa ($M_\oplus$)', fontsize=12)
    ax1.set_title(f'Crecimiento y Composición de un Embrión en {r_emb} AU\nDisco con Neptuno a 5.0 AU', fontsize=14)
    ax1.grid(True, which='both', ls='--', alpha=0.5)
    ax1.legend(fontsize=12)
    
    # ── Gráfico 2: Fracción de Hielo (%) ──
    ax2.plot(t_yr, frac_h2o, lw=3, color='dodgerblue')
    ax2.axhline(50, color='red', linestyle='--', alpha=0.5, label='Límite 50% (Puros Pebbles Helados)')
    ax2.axhline(10, color='green', linestyle='--', alpha=0.5, label='Límite 10% (Umbral Waterworld)')
    
    ax2.set_ylabel('Fracción de Hielo H$_2$O (%)', fontsize=12)
    ax2.set_xlabel('Tiempo (años)', fontsize=12)
    ax2.set_ylim(0, 60)
    ax2.grid(True, ls='--', alpha=0.5)
    ax2.legend(fontsize=12)
    
    plt.tight_layout()
    
    output_file = f"embryo_composition_{r_emb}au_nep.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Gráfico guardado exitosamente como '{output_file}' en la raíz del proyecto.")
    
    # Imprimir resumen final en consola
    print("\n" + "="*50)
    print("RESUMEN FINAL DEL EMBRIÓN")
    print("="*50)
    print(f"Masa Final Total : {M_total[-1]:.3f} M_earth")
    print(f"Masa Roca        : {M_silicates[-1]:.3f} M_earth")
    print(f"Masa Hielo H2O   : {M_h2o[-1]:.3f} M_earth")
    print(f"Fracción Hielo   : {frac_h2o[-1]:.1f} %")
    print("="*50)

if __name__ == "__main__":
    main()

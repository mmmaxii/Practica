import os
import numpy as np
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

def main():
    print("=" * 60)
    print("Ejemplo de uso: WaterworldPipeline con Init Automático")
    print("=" * 60)

    # 1. Creamos el pipeline y configuramos TODO en una sola línea
    # Le pasamos la POSICIÓN del gap (5.0 AU) para que refine la grilla ahí.
    pipeline = WaterworldPipeline(
        datadir="data/ejemplo_automatico",
        active_species=[],   # Especies químicas a inyectar
        grid_rmin=0.5 * c.au,            # Límite interno de la grilla
        grid_rmax=1000 * c.au,           # Límite externo
        Nr=500,                          # Número de celdas
        M_star_Msun=1.0,                 # 2 Masas Solares
        R_star_Rsun=1.0,                 # 2 Radios Solares
        T_star_K=4000.0,                # 10000 K
        gap_positions_au=[5.0],          # <- ¡AQUÍ ESTÁ! Le decimos dónde refinar
        snowline_au=2.0
    )

    # Nota: Hasta aquí el pipeline ya configuró la grilla, la estrella,
    # inicializó tripodpy, configuró la química, y las snowlines automáticamente.

    # 2. Agregamos la FISICA del gap (La profundidad / remoción de gas)
    # Como la grilla ya está fina en 5.0 AU, ahora le aplicamos el efecto del Júpiter.
    pipeline.setup_gap_duffell(
        M_planet=1.0 * 317.8 * c.M_earth, # Masa de 1 Júpiter
        a_planet_au=5.0
    )

    # 3. Actualizamos el simulador para que registre la física del gap
    pipeline.sim.update()

    # 4. Corremos la simulación
    print("\nLanzando integración temporal...")
    pipeline.run_integration(
        t_start_years=100,      # Comenzar snapshots en 100 años
        t_end_years=1e5,        # Terminar en 100,000 años
        num_snapshots=30        # Guardar 30 archivos HDF5
    )

    print("\n[OK] Simulación de ejemplo finalizada exitosamente.")

if __name__ == "__main__":
    main()

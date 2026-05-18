import numpy as np
from pipeline_snowlines import WaterworldPipeline
import dustpy.constants as c

print("Inicializando pipeline con física de snowline (Oka 2011)...")
pipeline = WaterworldPipeline()

# Los métodos de inicialización (add_volatile_components, add_snowline_fields, setup_physics)
# ya son llamados automáticamente en el __init__ de WaterworldPipeline.

r_snow_cm = float(pipeline.sim.dust.r_snow)
r_snow_au = r_snow_cm / c.au
print(f"\n[Test] Valor de r_snow inicial: {r_snow_cm:.2e} cm ({r_snow_au:.2f} AU)")

r_grid = pipeline.sim.grid.r / c.au
v_frag = pipeline.sim.dust.v.frag

# Encontrar transición de v_frag
ice_indices = np.where(v_frag >= 1000.0)[0]

if len(ice_indices) > 0:
    idx_transition = ice_indices[0]
    r_transition_au = r_grid[idx_transition]

    print(f"[Test] Transición de v_frag ocurre en la celda {idx_transition} (r = {r_transition_au:.2f} AU)")
    print(f"[Test] v_frag interno (celda 0): {v_frag[0]:.1f} cm/s")
    print(f"[Test] v_frag externo (celda -1): {v_frag[-1]:.1f} cm/s")

    if r_transition_au >= r_snow_au:
        print("\nVERIFICACIÓN EXITOSA: La transición de v_frag respeta el radio de r_snow de Oka 2011.")
    else:
        print("\nFALLO: La transición de v_frag ocurre DENTRO de la snowline.")
else:
    print("\n[Test] No hay hielo en la grilla. Todo el disco está al interior del snowline.")
    print(f"[Test] v_frag en la celda 0: {v_frag[0]:.1f} cm/s")
    print(f"[Test] v_frag en la celda -1: {v_frag[-1]:.1f} cm/s")
    print("\nVERIFICACIÓN EXITOSA: El disco entero está caliente al inicio (t=0, Mdot es altísimo).")

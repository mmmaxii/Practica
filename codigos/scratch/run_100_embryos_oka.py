import numpy as np
import dustpy.constants as c
import matplotlib.pyplot as plt
import sys
import os

# Asegurar que python encuentre los módulos en la carpeta raíz (PA3Py, etc)
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

datadir = "output_10kyr_oka"

print(f"Cargando simulaciones desde {datadir}...")
pa3 = PebbleAccretionModule3.from_datadir(datadir, M_star=1.0)

# 100 embriones desde 1 hasta 10 AU
embryo_locations_AU = np.linspace(0.5, 1.5, 50)
M0_g = 0.1 * c.M_earth  # 0.001 Masas terrestres

print(f"Plantando {len(embryo_locations_AU)} embriones con masa inicial {M0_g/c.M_earth} M_E...")
results = pa3.run_growth(embryo_locations_AU, M0_g=M0_g)

print("\n--- Resumen de Acreción (Todos los embriones) ---")
pa3.summary(results)

final_masses = []
for r in embryo_locations_AU:
    hist = results[r]
    if len(hist) > 0:
        final_masses.append(hist[-1][1] / c.M_earth)
    else:
        final_masses.append(M0_g / c.M_earth)

plt.figure(figsize=(8, 5))
plt.plot(embryo_locations_AU, final_masses, 'o-', color='teal', markersize=4, alpha=0.8)
plt.axhline(0.001, color='gray', ls='--', lw=1.5, label='Masa inicial (0.001 $M_\oplus$)')
plt.xlabel("Radio [AU]")
plt.ylabel("Masa Final [$M_\oplus$]")
plt.yscale("log")
plt.title("Crecimiento de 100 embriones - Snowline Dinámico Oka et al.")
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()

output_fig = "embryos_oka_1_10AU.png"
plt.savefig(output_fig, dpi=200)
print(f"Gráfico completo guardado en: {output_fig}")

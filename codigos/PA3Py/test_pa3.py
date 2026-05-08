import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
from PA3Py.plot_pa3 import PA3Diagnostics
import numpy as np
import matplotlib.pyplot as plt

DATADIR = r"C:\astro\Codigos practica + docs + papers\codigos\data_gaps_pipeline\t_5e6\sinusoidal_A10_fuerte"

# Grilla densa de embriones — concentrada en la zona de la snowline del H2O
EMBRYOS = sorted(list(np.logspace(np.log10(1.2), np.log10(30.0), num=300)))

print("Cargando datos del discos...")
pam3 = PebbleAccretionModule3.from_datadir(DATADIR, M_star=1.0)

print("\nEjecutando integracion de crecimiento...")
results = pam3.run_growth(EMBRYOS, M0_g=float(pam3.M_EARTH * 1e-3))
pam3.summary(results)

# ── Diagnosticos ───────────────────────────────────────────────────────────────
print("\nGenerando plots...")
diag = PA3Diagnostics(pam3, results, savedir="figs_pa3_test")

# Embriones representativos para plots temporales (composicion)
r_comp = [r for r in EMBRYOS if any(
    abs(r - x) < 0.15 for x in [1.5, 2.2, 3.5, 7.0]
)][:4]

diag.plot_waterworld_map()
diag.plot_with_disk_temperature()
diag.plot_composition(r_selected=r_comp)
diag.plot_hovmoller_mass()

plt.show()

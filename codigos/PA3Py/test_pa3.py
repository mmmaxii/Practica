import sys, os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
import numpy as np

DATADIR = r"C:\astro\Codigos practica + docs + papers\codigos\data_post_pipeline\pipeline_v3_Sigma_update"
EMBRYOS = [1.0, 3.0, 5.0, 10.0, 20.0, 30.0]

print("Cargando modulo 3...")
pam3 = PebbleAccretionModule3.from_datadir(DATADIR, M_star=1.0)
print("\nEjecutando integracion de crecimiento...")
results = pam3.run_growth(EMBRYOS, M0_g=float(pam3.M_EARTH * 1e-2)) # masa semilla ~ masa lunar / 100

pam3.summary(results)

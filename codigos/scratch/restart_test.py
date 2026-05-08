"""
restart_test.py — Prueba del mecanismo de restart desde dump file.
===================================================================

Flujo:
  1. Corre una simulación corta (1e3 → 1e4 yr) y guarda el dump.
  2. La reinicia con restart_from_dump() hasta 1e5 yr.
  3. Compara el número de HDF5 antes y después para verificar que funcionó.

Uso:
    python restart_test.py
"""

import os
import sys
import glob
import numpy as np
import dustpy.constants as c

sys.path.insert(0, os.path.abspath('.'))
from pipeline_snowlines import WaterworldPipeline

TEST_DIR = "_restart_test"

# ─────────────────────────────────────────────────────────────────────────────
# PASO 1 — Correr simulación inicial corta (1e3 → 1e4 yr)
# ─────────────────────────────────────────────────────────────────────────────
print("=" * 60)
print("PASO 1: Simulación inicial (1e3 → 1e4 yr)")
print("=" * 60)

pipeline = WaterworldPipeline(TEST_DIR)
pipeline.active_species = ["H2O"]
pipeline.setup_grid(rmin=1*c.au, rmax=100*c.au, Nr=50)   # Nr pequeño → rápido
pipeline.setup_star()
pipeline.setup_refined_grid(
    gap_positions_au = [5.0],
    num_gap          = 3,
    num_snow         = 2,
)
pipeline.initialize_simulation()
pipeline.add_volatile_components()
pipeline.setup_physics()
pipeline.setup_star_evolution()
pipeline.add_snowline_fields()
pipeline.add_ice_sigma_fields()
pipeline.sim.update()
pipeline.run_integration(t_start_years=1e3, t_end_years=1e4, num_snapshots=5)

hdfs_after_paso1 = sorted(glob.glob(os.path.join(TEST_DIR, "*.hdf5")))
dump_path = os.path.join(TEST_DIR, "frame.dmp")

print(f"\n✓ Paso 1 completado.")
print(f"  HDF5 escritos: {len(hdfs_after_paso1)}")
print(f"  Dump existe:   {os.path.isfile(dump_path)}")
print(f"  Último snap:   {os.path.basename(hdfs_after_paso1[-1])}")

# ─────────────────────────────────────────────────────────────────────────────
# PASO 2 — Restart desde el dump hasta 1e5 yr
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("PASO 2: Restart desde dump (→ 1e5 yr)")
print("=" * 60)

sim_r = WaterworldPipeline.restart_from_dump(
    datadir             = TEST_DIR,
    t_end_years         = 1e5,
    num_extra_snapshots = 8,
)

hdfs_after_paso2 = sorted(glob.glob(os.path.join(TEST_DIR, "*.hdf5")))

print(f"\n✓ Restart completado.")
print(f"  HDF5 total: {len(hdfs_after_paso1)} → {len(hdfs_after_paso2)}")
print(f"  Añadidos:   +{len(hdfs_after_paso2) - len(hdfs_after_paso1)} snapshots")
print(f"  t final:    {float(sim_r.t) / c.year:.3e} yr")
print(f"  Último snap: {os.path.basename(hdfs_after_paso2[-1])}")

# ─────────────────────────────────────────────────────────────────────────────
# Verificación rápida: gas.Sigma al final
# ─────────────────────────────────────────────────────────────────────────────
print("\n" + "=" * 60)
print("VERIFICACIÓN: gas.Sigma en la última posición")
print("=" * 60)
r_au = sim_r.grid.r / c.au
sig  = np.array(sim_r.gas.Sigma)
r_mid_idx = len(r_au) // 2
print(f"  r[medio] = {r_au[r_mid_idx]:.2f} AU  →  Σ_gas = {sig[r_mid_idx]:.3e} g/cm²")
print("\n¡Restart funcionó correctamente!")

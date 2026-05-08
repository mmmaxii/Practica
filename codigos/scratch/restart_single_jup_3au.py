"""
restart_single_jup_3au.py — Reinicio de single_jup_3.0au desde frame.dmp
=========================================================================

La simulación single_jup_3.0au quedó interrumpida con 55/60 snapshots.
Este script la reanuda desde el dump hasta completar 1e6 yr.

Uso:
    python restart_single_jup_3au.py
"""

import os
import glob
import sys
import dustpy.constants as c

sys.path.insert(0, os.path.abspath('.'))
from pipeline_snowlines import WaterworldPipeline

DATADIR  = os.path.join("data_1myr", "single_jup_3.0au")
T_END_YR = 1e6   # mismo t_end que el run original

# ── Diagnóstico previo ────────────────────────────────────────────────────────
hdfs_before = sorted(glob.glob(os.path.join(DATADIR, "*.hdf5")))
dump_path   = os.path.join(DATADIR, "frame.dmp")

print("=" * 60)
print("RESTART: single_jup_3.0au")
print("=" * 60)
print(f"  Directorio : {DATADIR}")
print(f"  HDF5 antes : {len(hdfs_before)} snapshots")
print(f"  Dump existe: {os.path.isfile(dump_path)}")
print(f"  t_end      : {T_END_YR:.1e} yr")
print()

if not os.path.isfile(dump_path):
    print("[ERROR] No se encontró el dump file. Abortando.")
    sys.exit(1)

# ── Restart ───────────────────────────────────────────────────────────────────
sim = WaterworldPipeline.restart_from_dump(
    datadir             = DATADIR,
    t_end_years         = T_END_YR,
    num_extra_snapshots = 30,   # por si el dump no cubre hasta 1e6 yr
)

# ── Reporte final ─────────────────────────────────────────────────────────────
hdfs_after = sorted(glob.glob(os.path.join(DATADIR, "*.hdf5")))

print()
print("=" * 60)
print("✓ Restart completado.")
print(f"  HDF5: {len(hdfs_before)} → {len(hdfs_after)}  (+{len(hdfs_after) - len(hdfs_before)} nuevos)")
print(f"  t final: {float(sim.t) / c.year:.3e} yr")
if hdfs_after:
    print(f"  Último snapshot: {os.path.basename(hdfs_after[-1])}")
print("=" * 60)

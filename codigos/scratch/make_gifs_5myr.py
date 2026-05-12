"""
make_gifs_5myr.py — Genera GIFs animados para todos los runs de data/5myr
==========================================================================

Para cada subdirectorio en data/5myr/ que tenga HDF5:
  - Crea data/5myr/<run>/gif/
  - Llama a disk_visualizer/disk_viz.py con los parámetros configurados
  - Guarda el GIF como data/5myr/<run>/gif/disco.gif

Uso
---
    py make_gifs_5myr.py
    py make_gifs_5myr.py --every 3 --fps 8 --rmax 25
    py make_gifs_5myr.py --runs single_jup_3.0au single_sat_5.0au
    py make_gifs_5myr.py --skip-existing

Argumentos (todos opcionales)
------------------------------
  --datadir     Carpeta raíz de datos      (default: data/5myr)
  --every       1 de cada N snapshots      (default: 2)
  --fps         Frames por segundo del GIF (default: 5)
  --rmax        Radio máximo del plot [AU] (default: 30)
  --dpi         DPI de los frames          (default: 100)
  --gif-name    Nombre del archivo GIF     (default: disco.gif)
  --runs        Lista de runs a procesar   (default: todos)
  --skip-existing  Omitir si el GIF ya existe
"""

import os
import sys
import glob
import argparse
import subprocess
import time

# ── Ruta al visualizador ──────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DISK_VIZ   = os.path.join(SCRIPT_DIR, "disk_visualizer", "disk_viz.py")


def find_runs(datadir):
    """Devuelve lista de subdirectorios de datadir que contengan HDF5."""
    runs = []
    for entry in sorted(os.scandir(datadir), key=lambda e: e.name):
        if not entry.is_dir():
            continue
        hdfs = glob.glob(os.path.join(entry.path, "data*.hdf5"))
        if hdfs:
            runs.append(entry.name)
        else:
            print(f"  [!] {entry.name}: sin HDF5 — omitido")
    return runs


def make_gif(run_name, datadir, gif_name, every, fps, rmax, dpi, skip_existing):
    """Genera el GIF para un run dado."""
    run_dir = os.path.join(datadir, run_name)
    gif_dir = os.path.join(run_dir, "gif")
    gif_path = os.path.join(gif_dir, gif_name)

    if skip_existing and os.path.isfile(gif_path):
        print(f"  [SKIP] {run_name}: GIF ya existe ({gif_path})")
        return True

    os.makedirs(gif_dir, exist_ok=True)

    cmd = [
        sys.executable, DISK_VIZ,
        "--datadir", run_dir,
        "--gif",     gif_path,
        "--every",   str(every),
        "--fps",     str(fps),
        "--rmax",    str(rmax),
        "--dpi",     str(dpi),
    ]

    print(f"\n{'='*70}")
    print(f"[GIF] {run_name}")
    print(f"  → Salida: {gif_path}")
    print(f"  → Comando: {' '.join(cmd)}")
    print(f"{'='*70}")

    t0 = time.time()
    result = subprocess.run(cmd, text=True)
    elapsed = (time.time() - t0) / 60

    if result.returncode == 0:
        size_mb = os.path.getsize(gif_path) / 1e6 if os.path.isfile(gif_path) else 0
        print(f"  ✓ Completado en {elapsed:.1f} min  |  {size_mb:.1f} MB")
        return True
    else:
        print(f"  ✗ ERROR en {run_name} (código {result.returncode})")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Genera GIFs para todos los runs de data/5myr."
    )
    parser.add_argument("--datadir",      default="data/5myr",
                        help="Carpeta raíz de datos (default: data/5myr)")
    parser.add_argument("--every",        type=int,   default=2,
                        help="1 de cada N snapshots (default: 2)")
    parser.add_argument("--fps",          type=int,   default=5,
                        help="FPS del GIF (default: 5)")
    parser.add_argument("--rmax",         type=float, default=30,
                        help="Radio máximo del plot en AU (default: 30)")
    parser.add_argument("--dpi",          type=int,   default=100,
                        help="DPI de los frames (default: 100)")
    parser.add_argument("--gif-name",     default="disco.gif",
                        help="Nombre del archivo GIF (default: disco.gif)")
    parser.add_argument("--runs",         nargs="+",  default=None,
                        help="Lista de runs a procesar (default: todos)")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Omitir runs cuyo GIF ya existe")
    args = parser.parse_args()

    # ── Validar paths ─────────────────────────────────────────────────────────
    if not os.path.isdir(args.datadir):
        print(f"[ERROR] Directorio no encontrado: {args.datadir}")
        sys.exit(1)
    if not os.path.isfile(DISK_VIZ):
        print(f"[ERROR] disk_viz.py no encontrado: {DISK_VIZ}")
        sys.exit(1)

    # ── Seleccionar runs ──────────────────────────────────────────────────────
    if args.runs:
        runs = args.runs
        # Verificar que existen
        for r in runs:
            if not os.path.isdir(os.path.join(args.datadir, r)):
                print(f"[ERROR] Run no encontrado: {os.path.join(args.datadir, r)}")
                sys.exit(1)
    else:
        print(f"\nBuscando runs en: {args.datadir}/")
        runs = find_runs(args.datadir)

    if not runs:
        print("[ERROR] No se encontraron runs con HDF5.")
        sys.exit(1)

    print(f"\n{'='*70}")
    print(f"  Runs a procesar : {len(runs)}")
    print(f"  every={args.every}  fps={args.fps}  rmax={args.rmax} AU  dpi={args.dpi}")
    print(f"  GIF → <run>/gif/{args.gif_name}")
    print(f"  skip-existing   : {args.skip_existing}")
    print(f"{'='*70}\n")

    # ── Procesar ──────────────────────────────────────────────────────────────
    ok, fail = [], []
    for i, run in enumerate(runs, 1):
        print(f"\n--- [{i}/{len(runs)}] {run} ---")
        success = make_gif(
            run_name      = run,
            datadir       = args.datadir,
            gif_name      = args.gif_name,
            every         = args.every,
            fps           = args.fps,
            rmax          = args.rmax,
            dpi           = args.dpi,
            skip_existing = args.skip_existing,
        )
        (ok if success else fail).append(run)

    # ── Resumen ───────────────────────────────────────────────────────────────
    print(f"\n{'='*70}")
    print(f"RESUMEN: {len(ok)}/{len(runs)} GIFs generados correctamente")
    if fail:
        print(f"  Fallidos ({len(fail)}):")
        for r in fail:
            print(f"    - {r}")
    print(f"{'='*70}\n")


if __name__ == "__main__":
    main()

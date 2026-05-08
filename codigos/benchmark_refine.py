"""
benchmark_refine.py — Visualiza la grilla antes y después de setup_refined_grid().
Patrón correcto (dustpylib): refinement ANTES de initialize_simulation().
"""
import sys, os
sys.path.insert(0, os.path.abspath('.'))

import numpy as np
import matplotlib.pyplot as plt
import dustpy.constants as c
from pipeline_snowlines import WaterworldPipeline

# ── Configuración de prueba ─────────────────────────────────────────────────
GAP_AU      = 5.0   # posición del gap
SNOWLINE_AU = 2.0   # posición de la snowline H2O
NR_BASE     = 100
RMIN        = 1.0 * c.au
RMAX        = 100.0 * c.au

# ── Pipeline: patrón correcto (dustpylib) ────────────────────────────────────
# ORDEN: setup_grid → setup_star → setup_refined_grid → initialize_simulation
pipeline = WaterworldPipeline("_bench_tmp")
pipeline.active_species = ["H2O"]
pipeline.setup_grid(rmin=RMIN, rmax=RMAX, Nr=NR_BASE)
pipeline.setup_star()

# ri ANTES: lo construimos igual que lo haría sim.initialize()
ri_before = np.geomspace(RMIN, RMAX, NR_BASE + 1)

# Refinar ANTES de initialize (el Field sim.grid.ri es constante después)
pipeline.setup_refined_grid(
    gap_positions_au = [GAP_AU],
    num_gap          = 3,
    num_snow         = 2,
    snowline_au      = SNOWLINE_AU,
)

pipeline.initialize_simulation()

# ri DESPUÉS: ahora sim.grid.ri tiene la grilla refinada
ri_after = np.array(pipeline.sim.grid.ri).copy()

print(f"\nNr antes:   {len(ri_before)-1}")
print(f"Nr después: {len(ri_after)-1}")
print(f"Celdas añadidas: {(len(ri_after)-len(ri_before))}")

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=False)
fig.suptitle("Benchmark setup_refined_grid() — Densidad de celdas radiales",
             fontsize=13, fontweight="bold")

r_before_au  = 0.5 * (ri_before[:-1] + ri_before[1:]) / c.au
dr_before_au = np.diff(ri_before) / c.au
r_after_au   = 0.5 * (ri_after[:-1]  + ri_after[1:])  / c.au
dr_after_au  = np.diff(ri_after)  / c.au

# ── Panel superior: dr(r) — toda la grilla ───────────────────────────────────
ax = axes[0]
ax.semilogy(r_before_au, dr_before_au, 'o-', ms=3, lw=1.2,
            color='steelblue', alpha=0.7, label=f"Antes  (Nr={len(r_before_au)})")
ax.semilogy(r_after_au,  dr_after_au,  'o-', ms=2, lw=1.2,
            color='tomato',    alpha=0.9, label=f"Después (Nr={len(r_after_au)})")
ax.axvline(SNOWLINE_AU, color='cyan',   lw=1.5, ls='--',
           label=f"H₂O snowline ({SNOWLINE_AU} AU)  num=2")
ax.axvline(GAP_AU,      color='orange', lw=1.5, ls='--',
           label=f"Gap ({GAP_AU} AU)  num=3")
ax.set_ylabel("Δr por celda [AU]")
ax.set_xlabel("Radio [AU]")
ax.set_title("Spacing radial entre celdas (valores menores = mayor resolución)")
ax.legend(fontsize=9)
ax.grid(True, which='both', alpha=0.3)

# ── Panel inferior: event plot zoom 1–15 AU ──────────────────────────────────
ax2 = axes[1]
ax2.eventplot(ri_before / c.au, lineoffsets=0, linelengths=0.7,
              color='steelblue', alpha=0.6, linewidths=0.9,
              label=f"Antes  ({len(ri_before)} interfaces)")
ax2.eventplot(ri_after  / c.au, lineoffsets=1, linelengths=0.7,
              color='tomato',    alpha=0.9, linewidths=0.9,
              label=f"Después ({len(ri_after)} interfaces)")
ax2.axvline(SNOWLINE_AU, color='cyan',   lw=1.5, ls='--')
ax2.axvline(GAP_AU,      color='orange', lw=1.5, ls='--')
ax2.set_xlim(0.8, 15.0)
ax2.set_yticks([0, 1])
ax2.set_yticklabels(["Antes", "Después"])
ax2.set_xlabel("Radio [AU]  (zoom 1–15 AU)")
ax2.set_title("Posición de interfaces (concentración visible cerca de 2 AU y 5 AU)")
ax2.legend(fontsize=9)
ax2.grid(True, which='both', alpha=0.3)

plt.tight_layout()
outfile = "figures/misc/benchmark_refine.png"
os.makedirs("figures/misc", exist_ok=True)
plt.savefig(outfile, dpi=150, bbox_inches='tight')
plt.show()
print(f"Plot guardado: {outfile}")

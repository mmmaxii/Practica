"""
test_pebble_accretion.py
========================
Prueba el módulo PebbleAccretion2.py con datos reales de tripodpy.

Uso:
    python test_pebble_accretion.py
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")          # sin ventana interactiva
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# ── Ajustar path ─────────────────────────────────────────────────────────
sys.path.insert(0, os.path.dirname(__file__))
from PebbleAccretion2 import PebbleAccretionModule

# ── Configuración ─────────────────────────────────────────────────────────
DATADIR   = r"C:\astro\Codigos practica + docs + papers\codigos\data_post_pipeline\pipeline_icefrac"
SAVEDIR   = r"C:\astro\Codigos practica + docs + papers\codigos\figs_pebbles"
M_STAR    = 1.0   # M☉

# Embriones a simular [AU] — desde interior al exterior, cruzando snowlines
EMBRYOS_AU = [1.0, 2.0, 3.0, 5.0, 8.0, 15.0, 30.0, 60.0]

# Masa inicial del embrión: ~1 M_luna ≈ 7.3e25 g  (planeta embrionario típico)
M0_g = 7.3e25

# ─────────────────────────────────────────────────────────────────────────
os.makedirs(SAVEDIR, exist_ok=True)

# ── 1. Cargar datos ───────────────────────────────────────────────────────
print("=" * 60)
print("  TEST: PebbleAccretionModule con pipeline_1000au")
print("=" * 60)

pam = PebbleAccretionModule.from_datadir(DATADIR, M_star=M_STAR, t_min_yr=0)

print(f"\nGrid radial: {pam.r[0]/pam.AU:.2f} – {pam.r[-1]/pam.AU:.2f} AU  (Nr={len(pam.r)})")
print(f"Tiempos:     {pam.times[0]/3.156e7:.1e} – {pam.times[-1]/3.156e7:.1e} yr  (Nt={pam.Nt})")

# Mostrar snowlines en el último snapshot
print("\nSnowlines en el último snapshot:")
for sp in ('H2O', 'CO2', 'CO'):
    r_s = pam.rsnow[sp][-1] / pam.AU
    print(f"  rsnow_{sp} = {r_s:.2f} AU")

# ── 2. Correr acreción ────────────────────────────────────────────────────
print(f"\nSimulando {len(EMBRYOS_AU)} embriones: {EMBRYOS_AU} AU")
results = pam.run_growth(EMBRYOS_AU, M0_g=M0_g)

# ── 3. Tabla resumen ──────────────────────────────────────────────────────
pam.summary(results)

# ── 4. Plot 1: M_core(t) para todos los embriones ─────────────────────────
fig, ax = plt.subplots(figsize=(10, 5))
colors = cm.plasma(np.linspace(0.1, 0.9, len(EMBRYOS_AU)))

for (r_au, hist), col in zip(results.items(), colors):
    if len(hist) == 0:
        continue
    t_kyr = hist[:, 0] / 3.156e10   # s → kyr
    M_Me  = hist[:, 1] / pam.M_EARTH
    ax.loglog(t_kyr, M_Me, lw=2, color=col, label=f"{r_au} AU")

# Línea de masa de aislamiento típica
ax.axhline(20 * (0.04/0.05)**3, ls=':', color='gray', lw=1.2, label=r'$M_{\rm iso}$ típica')
ax.set_xlabel("Tiempo [kyr]", fontsize=11)
ax.set_ylabel(r"$M_{\rm core}$ [$M_\oplus$]", fontsize=11)
ax.set_title("Crecimiento de núcleos por acreción de pebbles", fontsize=12)
ax.legend(fontsize=8, ncol=2, loc='upper left')
ax.grid(True, which='both', alpha=0.3)
fig.tight_layout()
path1 = os.path.join(SAVEDIR, "test_core_growth.pdf")
fig.savefig(path1, bbox_inches='tight')
print(f"\n→ Guardado: {path1}")

# ── 5. Plot 2: Composición final (stacked bar) ────────────────────────────
r_vals, f_H2O_arr, f_CO2_arr, f_sil_arr = [], [], [], []
for r_au, hist in results.items():
    if len(hist) == 0:
        continue
    _, M, H2O, CO2, sil, _, _ = hist[-1]
    r_vals.append(r_au)
    f_H2O_arr.append(100 * H2O / (M + 1e-30))
    f_CO2_arr.append(100 * CO2 / (M + 1e-30))
    f_sil_arr.append(100 * sil / (M + 1e-30))

fig2, ax2 = plt.subplots(figsize=(9, 5))
x = np.arange(len(r_vals))
w = 0.6
b1 = ax2.bar(x, f_sil_arr,  w, label="Silicatos",  color="#C0392B", alpha=0.9)
b2 = ax2.bar(x, f_CO2_arr,  w, bottom=f_sil_arr,
             label="CO₂ ice", color="#27AE60", alpha=0.9)
base3 = [a+b for a, b in zip(f_sil_arr, f_CO2_arr)]
b3 = ax2.bar(x, f_H2O_arr,  w, bottom=base3,
             label="H₂O ice", color="#2980B9", alpha=0.9)

# Línea de waterworld threshold
ax2.axhline(10, ls='--', color='k', lw=1.2, label='10% H₂O → Waterworld')

ax2.set_xticks(x)
ax2.set_xticklabels([f"{r} AU" for r in r_vals], rotation=30)
ax2.set_ylabel("Fracción en masa final [%]", fontsize=11)
ax2.set_title("Composición final del núcleo por posición del embrión", fontsize=12)
ax2.legend(fontsize=9, loc="upper right")
ax2.set_ylim(0, 110)
fig2.tight_layout()
path2 = os.path.join(SAVEDIR, "test_composition_bar.pdf")
fig2.savefig(path2, bbox_inches='tight')
print(f"→ Guardado: {path2}")

# ── 6. Plot 3: Mapa de crecimiento vs r y t ───────────────────────────────
fig3, axes3 = plt.subplots(1, 2, figsize=(13, 5))

for r_au, hist in results.items():
    if len(hist) < 2:
        continue
    t_kyr = hist[:, 0] / 3.156e10
    M_Me  = hist[:, 1] / pam.M_EARTH
    f_H2O = 100 * hist[:, 2] / (hist[:, 1] + 1e-30)

    c = cm.plasma(r_au / max(EMBRYOS_AU))
    axes3[0].loglog(t_kyr, M_Me, lw=2, color=c, label=f"{r_au} AU")
    axes3[1].semilogx(t_kyr, f_H2O, lw=2, color=c)

axes3[0].set_xlabel("t [kyr]"); axes3[0].set_ylabel(r"$M_{\rm core}$ [$M_\oplus$]")
axes3[0].set_title("Crecimiento del núcleo"); axes3[0].legend(fontsize=7, ncol=2)
axes3[0].grid(True, which='both', alpha=0.3)

axes3[1].axhline(10, ls='--', color='gray', lw=1.2, label='10% → Waterworld')
axes3[1].set_xlabel("t [kyr]"); axes3[1].set_ylabel("Fracción H₂O [%]")
axes3[1].set_title("Evolución de la fracción de agua")
axes3[1].legend(fontsize=8); axes3[1].grid(True, alpha=0.3)

# Colorbar manual por radio
sm = plt.cm.ScalarMappable(cmap='plasma',
     norm=plt.Normalize(vmin=min(EMBRYOS_AU), vmax=max(EMBRYOS_AU)))
fig3.colorbar(sm, ax=axes3, label="Radio del embrión [AU]", pad=0.02)
fig3.suptitle("Módulo de Acreción de Pebbles — pipeline_1000au", fontsize=13, y=1.01)
fig3.tight_layout()
path3 = os.path.join(SAVEDIR, "test_growth_composition.pdf")
fig3.savefig(path3, bbox_inches='tight')
print(f"→ Guardado: {path3}")

print("\n✅ Prueba completada.\n")

import sys, os
sys.path.insert(0, r'c:\astro\Codigos practica + docs + papers\codigos')
from PebbleAccretion2 import PebbleAccretionModule
import numpy as np

DATADIR = r"C:\astro\Codigos practica + docs + papers\codigos\data_post_pipeline\pipeline_v3_Sigma_update"
EMBRYOS = [1.0, 2.0, 3.0, 4.0, 5.0, 8.0, 12.0, 15.0, 20.0, 30.0, 50.0]

pam = PebbleAccretionModule.from_datadir(DATADIR, M_star=1.0)

# Verificar cuanta polvo hay en diferentes radios en t=0 y t=final
print("Verificando dust.Sigma[0, r, 1] (bin grande) en radios clave:")
for r_target_au in [1, 3, 5, 10, 20, 30, 50]:
    r_cm = r_target_au * pam.AU
    idx_r = np.argmin(np.abs(pam.r - r_cm))
    sig_t0  = pam.dust['Sigma'][0,  idx_r, 1]
    sig_tf  = pam.dust['Sigma'][-1, idx_r, 1]
    ice_t0  = pam.comp['H2O'][0,  idx_r] if pam._has_comp_sigma else -1
    ice_tf  = pam.comp['H2O'][-1, idx_r] if pam._has_comp_sigma else -1
    T_t0    = pam.gas['T'][0, idx_r]
    T_tf    = pam.gas['T'][-1, idx_r]
    print(f"  r={r_target_au:5.1f} AU | "
          f"Sigma_peb t0={sig_t0:.2e} tf={sig_tf:.2e} g/cm2 | "
          f"SigmaDust_H2O t0={ice_t0:.2e} tf={ice_tf:.2e} | "
          f"T t0={T_t0:.1f} tf={T_tf:.1f} K")

print()
results = pam.run_growth(EMBRYOS, M0_g=7.3e25)

out = []
out.append(f"{'r [AU]':>8} {'M_final':>10} {'H2O%':>8} {'CO2%':>8} {'Sil%':>8}  Tipo")
out.append("-" * 68)
for r_au, hist in results.items():
    if len(hist) == 0:
        out.append(f"{r_au:>8.1f}  sin accrecion")
        continue
    _, M, H2O, CO2, sil, rsnow, Miso = hist[-1]
    M_earth = M / pam.M_EARTH
    M0_earth = 7.3e25 / pam.M_EARTH
    fH2O = 100 * H2O / (M + 1e-30)
    fCO2 = 100 * CO2 / (M + 1e-30)
    fsil = 100 * sil / (M + 1e-30)
    accreted = M_earth - M0_earth
    tipo = "WATERWORLD" if fH2O > 10 else "Rocoso"
    out.append(f"{r_au:>8.1f}  {M_earth:>9.4f}  {fH2O:>7.2f}  {fCO2:>7.2f}  {fsil:>7.2f}  {tipo}  (acr={accreted:.4f} Me)")

out.append("")
out.append("Snowlines en t_final:")
for sp in ('H2O', 'CO2', 'CO'):
    r_s = pam.rsnow[sp][-1] / pam.AU
    out.append(f"  rsnow_{sp} = {r_s:.2f} AU")

out.append(f"Componentes: {'SigmaDust real (HDF5)' if pam._has_comp_sigma else 'Snowline fallback'}")
out.append(f"t range: {pam.times[0]/3.156e7:.1e} - {pam.times[-1]/3.156e7:.1e} yr")

result_text = "\n".join(out)
print(result_text)

# Guardar a archivo limpio
with open("pebble_results.txt", "w", encoding="ascii", errors="replace") as f:
    f.write(result_text)
print("\nGuardado: pebble_results.txt")

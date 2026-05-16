"""
PA3Py/PebbleAccretion3.py — Módulo de acreción de pebbles (Ormel 2017 & Drążkowska et al. 2023)
=====================================================================================

Física implementada y referencias:
  - M_PA_onset = St * η³ * M_star             (Drążkowska et al. 2023 Eq. 3, Ormel 2017 Eq 7.11)
  - M_hw/sh   = v_hw³ / (8 G Ω_K t_stop)      (Ormel 2017 Eq 7.9)
  - Ṁ_2D_hw   = √(8GM t_stop v_hw) Σ_peb      (Ormel 2017 Eq 7.13 Headwind)
  - Ṁ_2D_sh   = 2 R_H² Ω_K St^(2/3) Σ_peb     (Ormel 2017 Eq 7.13 Shear, Drążkowska et al. 2023 Eq 5)
  - Ṁ_3D      = 2π G M t_stop ρ_peb           (Ormel 2017 Eq 7.12)
  - Transición = Ṁ_2D * b_col / (b_col + H_peb √(8/π)) (Ormel 2017 Eq 7.24)
  - M_iso_peb = 25 M⊕ (H_gas/r / 0.05)³ (M_star/M_sun) (Drążkowska et al. 2023 Eq 6)
  - M < M_onset: Acreción de Safronov Balística (Ormel 2017 Eq 7.14)

Unidades internas: CGS (g, cm, s).
"""

import os
import glob
import h5py
import numpy as np
import dustpy.constants as c

class PebbleAccretionModule3:
    """
    Simula la acreción de pebbles sobre embriones planetarios usando
    los snapshots HDF5 de tripodpy. 
    Version actualizada a física pura Ormel 2017 y revisión PA.
    """

    # ── Constantes físicas (CGS) ──────────────────────────────────────────
    G_CGS   = c.G    # cm³ g⁻¹ s⁻²
    M_SUN   = c.M_sun    # g
    M_EARTH = c.M_earth    # g
    AU      = c.au # cm

    # Índice único de pebbles.
    # Todos los arrays de polvo (Nr, 5) deben compartir el mismo índice físico.
    # Los bins son: [a0, fudge*a1, a1, 0.5*amax, amax]. El último (-1) representa los pebbles.
    peb_idx = -1

    _ABUNDANCES = {'H2O': 1.6e-4, 'CO2': 4.0e-5, 'CO': 8.0e-5, 'silicates': 2.0e-3}

    # ══════════════════════════════════════════════════════════════════════
    # Constructor y carga HDF5
    # ══════════════════════════════════════════════════════════════════════

    @classmethod
    def from_datadir(cls, datadir, M_star=1.0, t_min_yr=0.0):
        files = sorted(glob.glob(os.path.join(datadir, 'data*.hdf5')))
        if not files:
            raise FileNotFoundError(f"No se encontraron archivos HDF5 en {datadir}")

        print(f"[from_datadir PA3Py] Leyendo {len(files)} snapshots desde {datadir}...")
        times_list, rsnow = [], {'H2O': [], 'CO2': [], 'CO': []}
        gas_S, gas_T, gas_cs, gas_eta, gas_nu, gas_alpha, gas_Hp = [], [], [], [], [], [], []
        dust_Sigma, dust_vr, dust_St, dust_H = [], [], [], []
        OmegaK_list = []   # leemos del HDF5 directamente
        comp_sigma = {}   
        r_grid = None

        for fpath in files:
            with h5py.File(fpath, 'r') as f:
                t_s = float(f['t'][()])
                if t_s < t_min_yr * 3.156e7:
                    continue

                times_list.append(t_s)

                if r_grid is None:
                    r_grid = f['grid/r'][:]
                    grid_keys = list(f['grid'].keys())
                    # Detectar TODAS las especies con SigmaDust_X, incluyendo silicates
                    detected  = [k.replace('SigmaDust_', '') for k in grid_keys if k.startswith('SigmaDust_')]
                    all_sps   = detected if detected else ['H2O', 'CO2', 'CO']
                    comp_sigma = {sp: [] for sp in all_sps}
                    print(f"[from_datadir PA3Py] Especies detectadas (volátiles + silicatos): {all_sps}")
                    print(f"[from_datadir PA3Py] Denominador de fracciones: dust.Sigma.sum(-1)")

                # Gas
                gas_S.append(f['gas/Sigma'][:])
                gas_T.append(f['gas/T'][:])
                gas_cs.append(f['gas/cs'][:])
                gas_eta.append(f['gas/eta'][:])
                gas_nu.append(f['gas/nu'][:])
                
                # Leemos alpha y Hp directamente en lugar de recalcularlos
                gas_alpha.append(f['gas/alpha'][:])
                gas_Hp.append(f['gas/Hp'][:])

                for sp in ('H2O', 'CO2', 'CO'):
                    key = f'grid/rsnow_{sp}'
                    rsnow[sp].append(float(f[key][()]) if key in f else np.nan)

                # Dust  — guardamos Sigma y otras propiedades (Nt, Nr, 5)
                dust_Sigma.append(f['dust/Sigma'][:])
                dust_vr.append(f['dust/v/rad'][:])
                dust_St.append(f['dust/St'][:])
                dust_H.append(f['dust/H'][:])

                # OmegaK — leer directamente del HDF5 (shape Nr,)
                if 'grid/OmegaK' in f:
                    OmegaK_list.append(f['grid/OmegaK'][:])

                # Química: todas las especies (volátiles + silicatos)
                for sp in comp_sigma:
                    new_key = f'grid/SigmaDust_{sp}'
                    arr = f[new_key][:] if new_key in f else np.zeros(len(r_grid))
                    comp_sigma[sp].append(arr)

        Nt = len(times_list)
        Nr = len(r_grid)
        print(f"[from_datadir PA3Py] {Nt} snapshots válidos  |  Nr = {Nr}")

        obj = cls.__new__(cls)
        obj.r = r_grid
        obj.times = np.array(times_list)
        obj.Nt = Nt
        obj.M_star_solar = M_star
        obj.M_star = M_star * cls.M_SUN

        for sp in ('H2O', 'CO2', 'CO'):
            rsnow[sp] = np.array(rsnow[sp])
        obj.rsnow = rsnow

        obj.gas = {
            'Sigma': np.array(gas_S),
            'T':     np.array(gas_T),
            'cs':    np.array(gas_cs),
            'eta':   np.array(gas_eta),
            'nu':    np.array(gas_nu),
            'alpha': np.array(gas_alpha),
            'Hp':    np.array(gas_Hp),
        }
        obj.dust = {
            'Sigma':       np.array(dust_Sigma),         # (Nt, Nr, 5)
            'Sigma_total': np.array(dust_Sigma).sum(-1), # (Nt, Nr)  ← denominador de fracciones
            'vr':          np.array(dust_vr),
            'St':          np.array(dust_St),
            'H':           np.array(dust_H),             # (Nt, Nr, 5)
        }

        if comp_sigma and all(len(v) > 0 for v in comp_sigma.values()):
            obj.comp = {sp: np.stack(comp_sigma[sp]) for sp in comp_sigma}  # (Nt, Nr) por especie
            obj._has_comp_sigma = True
            # Todas las especies detectadas (volátiles + silicatos)
            obj._volatile_sps = [sp for sp in comp_sigma]
            print(f"[from_datadir PA3Py] comp.keys = {list(obj.comp.keys())}")
        else:
            obj.comp = {}
            obj._has_comp_sigma = False
            obj._volatile_sps = []

        # OmegaK — usar datos del HDF5 si existen, si no calcular analíticamente
        if OmegaK_list:
            obj.Omega_K = np.array(OmegaK_list)   # (Nt, Nr) — del HDF5 (data.grid.OmegaK)
            obj._omegaK_2d = True
            print(f"[from_datadir PA3Py] OmegaK: leído del HDF5 — shape {obj.Omega_K.shape}")
        else:
            # Fallback analítico (1D, sin variación temporal)
            obj.Omega_K = np.sqrt(obj.G_CGS * obj.M_star / obj.r**3)  # (Nr,)
            obj._omegaK_2d = False
            print("[from_datadir PA3Py] OmegaK: calculado analíticamente (fallback)")
        
        # Asignamos H_gas y alpha directamente desde la lectura del HDF5
        # Evitamos cálculos manuales que puedan divergir de tripodpy
        obj.H_gas = obj.gas['Hp']
        obj.alpha = obj.gas['alpha']
        
        return obj

    def _comp_from_snowlines(self, t_idx, r_emb):
        fracs = {}
        for sp, abun in self._ABUNDANCES.items():
            if sp == 'CO':
                continue
            rsnow_t = self.rsnow.get(sp, np.array([np.nan]*self.Nt))
            r_snow  = float(rsnow_t[t_idx]) if not np.isnan(rsnow_t[t_idx]) else 0.0
            if sp == 'silicates' or r_emb > r_snow:
                fracs[sp] = abun
            else:
                fracs[sp] = 0.0
        total = sum(fracs.values())
        if total > 0:
            return {k: v / total for k, v in fracs.items()}
        return {'silicates': 1.0, 'H2O': 0.0, 'CO2': 0.0}

    # ══════════════════════════════════════════════════════════════════════
    # Helpers Interp
    # ══════════════════════════════════════════════════════════════════════

    def _interp(self, field_1d, r_emb):
        """
        Interpola logarítmicamente un campo radial 1D (definido en la grilla self.r)
        para obtener su valor en la posición específica del embrión.

        Parámetros:
        -----------
        field_1d : np.ndarray
            Arreglo 1D con los valores del campo en cada celda de la grilla radial (shape: Nr,).
        r_emb : float
            Posición radial del embrión en cm.

        Retorna:
        --------
        float
            Valor interpolado del campo en la posición r_emb.
        """
        return float(np.interp(np.log(r_emb), np.log(self.r), field_1d))

    def _local(self, t, r_emb):
        """
        Extrae e interpola las propiedades locales del disco (gas y polvo)
        en un instante de tiempo temporal y en la posición del embrión,
        que son necesarias para evaluar el flujo de acreción.

        Parámetros:
        -----------
        t : int
            Índice temporal correspondiente al snapshot actual (de 0 a Nt-1).
        r_emb : float
            Posición radial del embrión en cm.

        Retorna:
        --------
        dict
            Diccionario con las propiedades locales:
            - Sigma_peb (float): Densidad superficial de pebbles [g/cm²].
            - eta (float): Parámetro adimensional del gradiente de presión del gas.
            - alpha (float): Parámetro de viscosidad turbulenta de Shakura-Sunyaev.
            - St (float): Número de Stokes máximo del polvo (pebbles).
            - H_gas (float): Escala de altura del gas [cm].
            - H_peb (float): Escala de altura de los pebbles [cm].
            - Omega (float): Frecuencia Kepleriana local [rad/s].
            - v_K (float): Velocidad Kepleriana local [cm/s].
            - v_hw (float): Velocidad del headwind (viento de cara) de los pebbles [cm/s].
        """
        I = lambda arr: self._interp(arr[t], r_emb)
        
        # Usamos peb_idx de forma consistente para todas las propiedades (Nr, 5)
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, self.peb_idx], r_emb)
        eta   = I(self.gas['eta'])
        St    = self._interp(self.dust['St'][t, :, self.peb_idx], r_emb)
        H_gas = self._interp(self.H_gas[t], r_emb)
        
        # dust.H ya incluye sedimentación y mezcla turbulenta usando dust.delta.vert
        # Evitamos inconsistencias entre viscosidad (gas.alpha) y mezcla vertical.
        H_peb = self._interp(self.dust['H'][t, :, self.peb_idx], r_emb)
        
        # OmegaK: (Nt, Nr) si viene del HDF5, (Nr,) si es fallback analítico
        omega_1d = self.Omega_K[t] if self._omegaK_2d else self.Omega_K
        Omega = float(np.interp(np.log(r_emb), np.log(self.r), omega_1d))
        v_K   = Omega * r_emb
        v_hw  = eta * v_K
        return dict(Sigma_peb=Sigma_peb, eta=eta, alpha=I(self.alpha),
                    St=St, H_gas=H_gas, H_peb=H_peb, Omega=Omega, v_K=v_K, v_hw=v_hw)

    # ══════════════════════════════════════════════════════════════════════
    # Física Ormel 2017 & Drążkowska et al. 2023
    # ══════════════════════════════════════════════════════════════════════

    def _pebble_flux(self, t, r_emb):
        """
        Ṁ_peb = 2π r Σ_peb |v_r|
        Mantenemos lectura directa de tripodpy para conservación real de masa.
        """
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, self.peb_idx], r_emb)
        v_r       = self._interp(self.dust['vr'][t, :, self.peb_idx],   r_emb)
        return 2 * np.pi * r_emb * Sigma_peb * abs(v_r)

    def _isolation_mass(self, r_emb, t):
        """
        Aislamiento por Pebbles (Drążkowska et al. 2023 Eq. 6):
          M_iso_peb = 25 M_earth * (H_gas/r / 0.05)^3 * (M_star/M_sun)
        (Equivalente a Ormel 2017 Eq 7.17, con factor 25 ajustado)
        """
        p = self._local(t, r_emb)
        h_gas = p['H_gas'] / r_emb
        M_iso = 25 * (h_gas / 0.05)**3 * self.M_star_solar * self.M_EARTH
        return max(M_iso, 0.1 * self.M_EARTH)

    def _accretion_rate(self, M_core, r_emb, t):
        if M_core <= 0: return 0.0
        p = self._local(t, r_emb)
        if p['Sigma_peb'] < 1e-30: return 0.0

        G, M, Omega = self.G_CGS, M_core, p['Omega']
        St, v_hw, Sigma = p['St'], p['v_hw'], p['Sigma_peb']
        
        # Tiempo de frenado dimensional
        t_stop = St / Omega

        # Onset of Pebble Accretion (Eq. 3 / Eq. 7.11 sin 1/8)
        M_PA_onset = St * (p['eta']**3) * self.M_star

        # Capa de pebbles
        # Ya no calculamos H_peb manualmente. Se lee dust.H directamente que incluye la mezcla vertical.
        H_peb = max(p['H_peb'], 1e-10 * p['H_gas'])
        rho_peb = Sigma / (np.sqrt(2 * np.pi) * H_peb)

        # Régimen Pre-Pebble / Safronov Balístico (Ormel Ec 7.14)
        # DESACTIVADO por solicitud del usuario para forzar crecimiento rápido
        # if M < M_PA_onset:
        #     rho_core = 3.0  # g/cm3 asumiendo rocoso
        #     R_pl = (3 * M / (4 * np.pi * rho_core))**(1/3)
        #     # Factor de corrección a evitar divergencias extremas
        #     v_impact = max(v_hw, 1.0) 
        #     Mdot_Saf = (2 * np.pi * R_pl * G * M / v_impact) * rho_peb
        #     return Mdot_Saf

        # Masa de transición Headwind/Shear (Ormel Ec 7.9)
        M_hw_sh = (v_hw**3) / (8 * G * Omega * St)

        # Regímenes 2D (Ormel Ec 7.13)
        if M < M_hw_sh:
            Mdot_2D = np.sqrt(8 * G * M * t_stop * v_hw) * Sigma
            b_col   = np.sqrt(2 * G * M * t_stop / max(v_hw, 1e-5))  # b_hw
        else:
            R_H     = r_emb * (M / (3 * self.M_star))**(1/3)
            Mdot_2D = 2 * R_H**2 * Omega * St**(2/3) * Sigma
            b_col   = (St**(1/3)) * R_H                              # b_sh
            
        # Transición 2D-3D turbulencia suave (Ormel 7.24)
        denominator = b_col + H_peb * np.sqrt(8.0 / np.pi)
        Mdot_eff = Mdot_2D * (b_col / denominator)
        
        return max(Mdot_eff, 0.0)

    # ══════════════════════════════════════════════════════════════════════
    # API / Loop Principal
    # ══════════════════════════════════════════════════════════════════════

    def run_growth(self, embryo_locations_AU, M0_g=None):
        if M0_g is None:
            M0_g = 1e-3 * self.M_EARTH
        locs_outer_to_inner = sorted(embryo_locations_AU, reverse=True)
        M_core   = {r: float(M0_g)              for r in locs_outer_to_inner}
        M_X      = {r: {'H2O': 0.0, 'silicates': 0.0} for r in locs_outer_to_inner}
        active   = {r: True                      for r in locs_outer_to_inner}
        histories= {r: []                        for r in locs_outer_to_inner}

        for i in range(self.Nt):
            dt = float(self.times[i] - (self.times[i-1] if i > 0 else 0.0))
            if dt <= 0: continue

            r_snow_AU = float(self.rsnow['H2O'][i]) / self.AU \
                        if not np.isnan(self.rsnow.get('H2O', [np.nan])[i]) else np.nan

            flux_consumed = 0.0

            for r_au in locs_outer_to_inner:
                if not active[r_au]: continue
                r_emb = r_au * self.AU

                M_iso = self._isolation_mass(r_emb, i)
                if M_core[r_au] >= M_iso:
                    active[r_au] = False
                    continue

                Mdot_peb_disk  = self._pebble_flux(i, r_emb)
                Mdot_peb_avail = max(0.0, Mdot_peb_disk - flux_consumed)

                Mdot_core_r = self._accretion_rate(M_core[r_au], r_emb, i)
                Mdot_core_r = min(Mdot_core_r, Mdot_peb_avail)

                dM = Mdot_core_r * dt
                dM = min(dM, max(0.0, M_iso - M_core[r_au]))
                flux_consumed += Mdot_core_r

                if r_emb > (r_snow_AU * self.AU if not np.isnan(r_snow_AU) else 0.0):
                    # Fuera de la snowline: 50% H2O, 50% Silicatos
                    M_X[r_au]['H2O'] += 0.5 * dM
                    M_X[r_au]['silicates'] += 0.5 * dM
                else:
                    # Dentro de la snowline: 100% Silicatos
                    M_X[r_au]['silicates'] += dM

                M_core[r_au] += dM

                histories[r_au].append((
                    self.times[i], M_core[r_au], M_X[r_au].get('H2O', 0),
                    M_X[r_au].get('CO2', 0), M_X[r_au].get('silicates', 0),
                    r_snow_AU, M_iso,
                ))

        results = {
            r_au: (np.array(histories[r_au]) if histories[r_au] else np.empty((0, 7)))
            for r_au in embryo_locations_AU
        }
        return results

    def summary(self, results):
        """Imprime tabla resumen de composicion final con M_iso."""
        SEP = '-' * 80
        vol_sps = ['H2O']
        header_sps = '  '.join(f'{sp+"%":>7}' for sp in vol_sps)
        print(f"\n{SEP}")
        print(f"{'r [AU]':>8} {'M [ME]':>9} {'M_iso [ME]':>11}  {header_sps}  {'Sil%':>7}  Tipo")
        print(SEP)
        
        for r_au, hist in results.items():
            if len(hist) == 0:
                print(f"{r_au:>8.2f}  -- sin acrecion")
                continue
            _, M, *_, M_iso = hist[-1]
            M_total = M
            row_vals = {}
            for j, sp in enumerate(vol_sps):
                # columna j+2 en historial: H2O=2, CO2=3, ...
                col = 2 + j
                M_sp = hist[-1][col] if col < len(hist[-1]) else 0.0
                row_vals[sp] = 100 * M_sp / (M_total + 1e-30)
            # Silicatos = masa no asignada a volátiles
            M_vol = sum(hist[-1][2 + j] for j in range(len(vol_sps)) if 2 + j < len(hist[-1]))
            f_sil  = 100 * max(0.0, M_total - M_vol) / (M_total + 1e-30)
            f_h2o  = row_vals.get('H2O', 0.0)
            tipo   = "[W] Waterworld" if f_h2o > 10 else "[R] Rocoso"
            vals_str = '  '.join(f"{row_vals.get(sp, 0.0):>6.1f}" for sp in vol_sps)
            print(f"{r_au:>8.2f}  {M/self.M_EARTH:>8.3f}  "
                  f"{M_iso/self.M_EARTH:>10.2f}  {vals_str}  {f_sil:>6.1f}  {tipo}")
        print(f"{SEP}\n")





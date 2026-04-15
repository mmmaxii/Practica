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

import h5py
import numpy as np

class PebbleAccretionModule3:
    """
    Simula la acreción de pebbles sobre embriones planetarios usando
    los snapshots HDF5 de tripodpy. 
    Version actualizada a física pura Ormel 2017 y revisión PA.
    """

    # ── Constantes físicas (CGS) ──────────────────────────────────────────
    G_CGS   = 6.6743e-8    # cm³ g⁻¹ s⁻²
    M_SUN   = 1.9884e33    # g
    M_EARTH = 5.9722e27    # g
    AU      = 1.4959787e13 # cm

    _ABUNDANCES = {'H2O': 1.6e-4, 'CO2': 4.0e-5, 'CO': 8.0e-5, 'silicates': 2.0e-3}

    # ══════════════════════════════════════════════════════════════════════
    # Constructor y carga HDF5
    # ══════════════════════════════════════════════════════════════════════

    @classmethod
    def from_datadir(cls, datadir, M_star=1.0, t_min_yr=0.0):
        import glob, os

        files = sorted(glob.glob(os.path.join(datadir, 'data*.hdf5')))
        if not files:
            raise FileNotFoundError(f"No se encontraron archivos HDF5 en {datadir}")

        print(f"[from_datadir PA3Py] Leyendo {len(files)} snapshots desde {datadir}...")
        times_list, rsnow = [], {'H2O': [], 'CO2': [], 'CO': []}
        gas_S, gas_T, gas_cs, gas_eta, gas_nu = [], [], [], [], []
        dust_Sigma, dust_vr, dust_St = [], [], []
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
                    detected  = [k.replace('SigmaDust_', '') for k in grid_keys if k.startswith('SigmaDust_')]
                    vol_sps   = detected if detected else ['H2O', 'CO2', 'CO']
                    comp_sigma = {sp: [] for sp in vol_sps}
                    comp_sigma['silicates'] = []
                    print(f"[from_datadir PA3Py] Especies detectadas: {vol_sps} + silicates")

                # Gas
                gas_S.append(f['gas/Sigma'][:])
                gas_T.append(f['gas/T'][:])
                gas_cs.append(f['gas/cs'][:])
                gas_eta.append(f['gas/eta'][:])
                gas_nu.append(f['gas/nu'][:])

                for sp in ('H2O', 'CO2', 'CO'):
                    key = f'grid/rsnow_{sp}'
                    rsnow[sp].append(float(f[key][()]) if key in f else np.nan)

                # Dust
                dust_Sigma.append(f['dust/Sigma'][:])
                dust_vr.append(f['dust/v/rad'][:])
                dust_St.append(f['dust/St'][:])

                # Química
                dust_snap = f['dust/Sigma'][:].sum(-1)
                vol_total = np.zeros_like(dust_snap)

                for sp in [k for k in comp_sigma if k != 'silicates']:
                    new_key = f'grid/SigmaDust_{sp}'
                    if new_key in f:
                        arr = f[new_key][:]
                    else:
                        arr = np.zeros(len(r_grid))
                    comp_sigma[sp].append(arr)
                    vol_total += arr

                comp_sigma['silicates'].append(np.maximum(0.0, dust_snap - vol_total))

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
        }
        obj.dust = {
            'Sigma': np.array(dust_Sigma),
            'vr':    np.array(dust_vr),
            'St':    np.array(dust_St),
        }

        if comp_sigma and all(len(v) > 0 for v in comp_sigma.values()):
            obj.comp = {sp: np.stack(comp_sigma[sp]) for sp in comp_sigma}
            obj._has_comp_sigma = True
        else:
            obj.comp = {}
            obj._has_comp_sigma = False

        obj.Omega_K = np.sqrt(obj.G_CGS * obj.M_star / obj.r**3)
        obj.H_gas   = obj.gas['cs'] / obj.Omega_K[np.newaxis, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            alpha_raw = obj.gas['nu'] / (obj.gas['cs'] * obj.H_gas)
        obj.alpha = np.clip(np.where(np.isfinite(alpha_raw), alpha_raw, 1e-3), 1e-5, 1e-1)
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
        return float(np.interp(np.log(r_emb), np.log(self.r), field_1d))

    def _local(self, t, r_emb):
        I = lambda arr: self._interp(arr[t], r_emb)
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, 1], r_emb)
        eta   = I(self.gas['eta'])
        St    = self._interp(self.dust['St'][t, :, -1], r_emb)
        H_gas = self._interp(self.H_gas[t], r_emb)
        Omega = float(np.interp(np.log(r_emb), np.log(self.r), self.Omega_K))
        v_K   = Omega * r_emb
        v_hw  = eta * v_K
        return dict(Sigma_peb=Sigma_peb, eta=eta, alpha=I(self.alpha),
                    St=St, H_gas=H_gas, Omega=Omega, v_K=v_K, v_hw=v_hw)

    # ══════════════════════════════════════════════════════════════════════
    # Física Ormel 2017 & Drążkowska et al. 2023
    # ══════════════════════════════════════════════════════════════════════

    def _pebble_flux(self, t, r_emb):
        """
        Ṁ_peb = 2π r Σ_peb |v_r|
        Mantenemos lectura directa de tripodpy para conservación real de masa.
        """
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, 1], r_emb)
        v_r       = self._interp(self.dust['vr'][t, :, -1],   r_emb)
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

        # Capa de pebbles (Ec 7.25 Ormel)
        H_peb = p['H_gas'] * np.sqrt(p['alpha'] / (p['alpha'] + St))
        H_peb = max(H_peb, 1e-10 * p['H_gas'])
        rho_peb = Sigma / (np.sqrt(2 * np.pi) * H_peb)

        # Régimen Pre-Pebble / Safronov Balístico (Ormel Ec 7.14)
        if M < M_PA_onset:
            rho_core = 3.0  # g/cm3 asumiendo rocoso
            R_pl = (3 * M / (4 * np.pi * rho_core))**(1/3)
            # Factor de corrección a evitar divergencias extremas
            v_impact = max(v_hw, 1.0) 
            Mdot_Saf = (2 * np.pi * R_pl * G * M / v_impact) * rho_peb
            return Mdot_Saf

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

    def run_growth(self, embryo_locations_AU, M0_g=1e24):
        locs_outer_to_inner = sorted(embryo_locations_AU, reverse=True)
        M_core   = {r: float(M0_g)              for r in locs_outer_to_inner}
        M_X      = {r: {sp: 0.0 for sp in self.comp} for r in locs_outer_to_inner}
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

                if self._has_comp_sigma:
                    sig_ice_total = sum(self._interp(self.comp[sp][i], r_emb) for sp in self.comp)
                    if sig_ice_total > 1e-100:
                        for sp in self.comp:
                            frac = self._interp(self.comp[sp][i], r_emb)
                            frac = np.clip(frac / sig_ice_total, 0.0, 1.0)
                            M_X[r_au][sp] += frac * dM
                    else:
                        fracs = self._comp_from_snowlines(i, r_emb)
                        for sp, frac in fracs.items():
                            M_X[r_au][sp] = M_X[r_au].get(sp, 0.0) + frac * dM
                else:
                    fracs = self._comp_from_snowlines(i, r_emb)
                    for sp, frac in fracs.items():
                        M_X[r_au][sp] = M_X[r_au].get(sp, 0.0) + frac * dM

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
        print(f"\n{SEP}")
        print(f"{'r [AU]':>8} {'M [ME]':>9} {'M_iso [ME]':>11} "
              f"{'H2O%':>7} {'CO2%':>7} {'Sil%':>7}  Tipo")
        print(SEP)
        for r_au, hist in results.items():
            if len(hist) == 0:
                print(f"{r_au:>8.2f}  -- sin acrecion")
                continue
            _, M, H2O, CO2, sil, _, M_iso = hist[-1]
            f_H2O = 100 * H2O / (M + 1e-30)
            f_CO2 = 100 * CO2 / (M + 1e-30)
            f_sil = 100 * sil / (M + 1e-30)
            tipo  = "[W] Waterworld" if f_H2O > 10 else "[R] Rocoso"
            print(f"{r_au:>8.2f}  {M/self.M_EARTH:>8.3f}  "
                  f"{M_iso/self.M_EARTH:>10.2f}  {f_H2O:>6.1f}  "
                  f"{f_CO2:>6.1f}  {f_sil:>6.1f}  {tipo}")
        print(f"{SEP}\n")




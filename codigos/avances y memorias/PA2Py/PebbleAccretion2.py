"""
PebbleAccretion2.py — Módulo de acreción de pebbles post-evolución de disco
============================================================================

Física implementada y referencias:
  - delta_v = η v_K               (Lambrechts & Johansen 2012, A&A 544, A32, Eq. 16)
  - R_acc = min(R_H, b_Bondi, b_Hill):
      b_Bondi = √St · GM/Δv²      (Ormel & Klahr 2010, A&A 520, A43, Eq. A9)
      b_Hill  = (St/0.1)^(1/3) RH  (Lambrechts & Johansen 2012, Eq. A10)
  - H_peb = H_gas √(α/(α+St))    (Youdin & Lithwick 2007, Icarus 192, 588, Eq. 12)
  - Ṁ_2D  = 2 R_acc Σ_peb Δv    (Ormel & Klahr 2010)
  - Ṁ_3D  = π R_acc² ρ_peb Δv   (Ormel & Klahr 2010)
  - Transición 2D/3D via erf(z),  z = R_acc/(√2 H_peb)
                                   (Liu & Ormel 2018, A&A 615, A178, Appendix A.1)
  - M_iso = 25 M⊕(h/0.05)³ f_α f_P (Bitsch et al. 2018, A&A 612, A30, Eq. 5)
      f_α = 0.34(log10(α/1e-3)+3)+0.66
      f_P = 1 - (d lnP/d lnr + 2.5)/6,  d lnP/d lnr = -2η/h²

Unidades internas: CGS (g, cm, s).
"""

import h5py
import numpy as np
from scipy.special import erf


class PebbleAccretionModule:
    """
    Simula la acreción de pebbles sobre embriones planetarios usando
    los snapshots HDF5 de tripodpy. El disco se trata como fondo fijo
    en cada snapshot (válido para t_accretion << t_disk).
    """

    # ── Constantes físicas (CGS) ──────────────────────────────────────────
    G_CGS   = 6.674e-8     # cm³ g⁻¹ s⁻²
    M_SUN   = 1.989e33     # g
    M_EARTH = 5.972e27     # g
    AU      = 1.4959787e13 # cm

    # Abundancias de masa relativa al polvo sólido (solo fallback: se usa cuando
    # SigmaIce_X no está disponible en el HDF5).
    #
    # H2O corregido de 1.6e-4 → 9.0e-4 para ser consistente con pipeline_snowlines.py
    # (razón hielo:roca ~ 1:1, Drążkowska & Alibert 2017). Con 1.6e-4 la fracción
    # de agua resultante era ~7%, imposibilitando la clasificación de waterworld.
    _ABUNDANCES = {'H2O':1.6e-4, 'CO2': 4.0e-5, 'CO': 8.0e-5, 'silicates': 2.0e-3}

    # ══════════════════════════════════════════════════════════════════════
    # Constructor principal (desde un único HDF5 combinado)
    # ══════════════════════════════════════════════════════════════════════

    def __init__(self, h5file, M_star=1.0):
        """
        Parameters
        ----------
        h5file : str   — HDF5 generado por tripodpy (read.all()).
        M_star : float — Masa estelar [M☉]. Default: 1.0.
        """
        f = h5py.File(h5file, 'r')

        # ── Grid radial: manejar shape (Nt, Nr) o (Nr,) ─────────────────
        r_raw  = f['grid/r'][:]
        self.r = r_raw[0] if r_raw.ndim == 2 else r_raw  # (Nr,) [cm]

        # ── Tiempo y snowlines ────────────────────────────────────────────
        self.times = f['t'][:]                            # (Nt,) [s]
        self.Nt    = len(self.times)
        self.rsnow = {}
        for sp in ('H2O', 'CO2', 'CO'):
            key = f'grid/rsnow_{sp}'
            self.rsnow[sp] = f[key][:] if key in f else np.full(self.Nt, np.nan)

        # ── Gas — shape (Nt, Nr) ─────────────────────────────────────────
        self.gas = {k: f[f'gas/{k}'][:] for k in ('Sigma', 'T', 'cs', 'eta', 'nu')}

        # ── Polvo — Sigma(Nt,Nr,2), vr(Nt,Nr,5), St(Nt,Nr,5) ───────────
        self.dust = {
            'Sigma': f['dust/Sigma'][:],
            'vr':    f['dust/v.rad'][:],
            'St':    f['dust/St'][:],
        }

        # ── Componentes químicos — shape (Nt, Nr) ──────────────────────
        # Prioridad de lectura:
        #   1. grid/SigmaDust_{sp}     — campo real del componente (nuevo, add_ice_sigma_fields)
        #   2. grid/SigmaIce_{sp}      — campo estimado por fracción+theta (prev. versión)
        #   3. components/{sp}/dust/Sigma — formato original de tripodpy
        # Silicatos: calculados como polvo_total − Σ(volátiles), no leidos directamente.
        self.comp = {}
        grid_keys  = list(f['grid'].keys()) if 'grid' in f else []
        volatile_sps = [k.replace('SigmaDust_', '') for k in grid_keys
                        if k.startswith('SigmaDust_')]
        if not volatile_sps:  # fallback: lista estándar
            volatile_sps = ['H2O', 'CO2', 'CO']

        dust_total_snap = self.dust['Sigma'].sum(-1)   # (Nt, Nr)
        volatile_total  = np.zeros_like(dust_total_snap)

        for sp in volatile_sps:
            new_key = f'grid/SigmaDust_{sp}'
            ice_key = f'grid/SigmaIce_{sp}'
            old_key = f'components/{sp}/dust/Sigma'
            if new_key in f:
                arr = f[new_key][:]
            elif ice_key in f:
                arr = f[ice_key][:]
            elif old_key in f:
                arr = f[old_key][:]
            else:
                arr = np.zeros((self.Nt, len(self.r)))
            self.comp[sp] = arr
            volatile_total += arr

        # Silicatos como fracción remanente del polvo total
        self.comp['silicates'] = np.maximum(0.0, dust_total_snap - volatile_total)
        f.close()

        # ── Parámetros estelares y derivados del disco ───────────────────
        self.M_star  = M_star * self.M_SUN
        self.Omega_K = np.sqrt(self.G_CGS * self.M_star / self.r**3)  # (Nr,)

        # H_gas(t,r) = cs/Omega_K   →  shape (Nt, Nr)
        self.H_gas = self.gas['cs'] / self.Omega_K[np.newaxis, :]

        # α(t,r) = ν/(cs·H_gas)   →  shape (Nt, Nr)  [Youdin & Lithwick 2007]
        with np.errstate(divide='ignore', invalid='ignore'):
            alpha_raw = self.gas['nu'] / (self.gas['cs'] * self.H_gas)
        self.alpha = np.clip(np.where(np.isfinite(alpha_raw), alpha_raw, 1e-3),
                             1e-5, 1e-1)
        # Marcador: componentes no disponibles (se completarán en from_datadir)
        if not hasattr(self, 'comp'):
            self.comp = {}
        self._has_comp_sigma = bool(self.comp)

    # ══════════════════════════════════════════════════════════════════════
    # Constructor alternativo: apila múltiples HDF5 de tripodpy
    # ══════════════════════════════════════════════════════════════════════

    @classmethod
    def from_datadir(cls, datadir, M_star=1.0, t_min_yr=0.0):
        """
        Carga todos los snapshots HDF5 de un directorio de tripodpy y
        construye el módulo apilando los datos en arrays (Nt, Nr, ...).

        Parameters
        ----------
        datadir  : str   — carpeta con data0000.hdf5, data0001.hdf5, ...
        M_star   : float — masa estelar [M☉].
        t_min_yr : float — excluir snapshots anteriores a t_min_yr años.

        Returns
        -------
        PebbleAccretionModule listo para usar.
        """
        import glob, os

        files = sorted(glob.glob(os.path.join(datadir, 'data*.hdf5')))
        if not files:
            raise FileNotFoundError(f"No se encontraron archivos HDF5 en {datadir}")

        print(f"[from_datadir] Leyendo {len(files)} snapshots desde {datadir}...")

        # ── Arrays de acumulación ────────────────────────────────────────
        times_list, rsnow = [], {'H2O': [], 'CO2': [], 'CO': []}
        gas_S, gas_T, gas_cs, gas_eta, gas_nu = [], [], [], [], []
        dust_Sigma, dust_vr, dust_St = [], [], []
        comp_sigma = {}   # poblado dinámicamente tras leer el primer archivo
        r_grid = None

        for fpath in files:
            with h5py.File(fpath, 'r') as f:
                t_s = float(f['t'][()])
                if t_s < t_min_yr * 3.156e7:
                    continue

                times_list.append(t_s)

                if r_grid is None:
                    r_grid = f['grid/r'][:]   # (Nr,)
                    # ── Detectar especies volátiles en el primer archivo ──────
                    grid_keys = list(f['grid'].keys())
                    detected  = [k.replace('SigmaDust_', '') for k in grid_keys
                                 if k.startswith('SigmaDust_')]
                    vol_sps   = detected if detected else ['H2O', 'CO2', 'CO']
                    comp_sigma = {sp: [] for sp in vol_sps}
                    comp_sigma['silicates'] = []
                    print(f"[from_datadir] Especies detectadas: {vol_sps} + silicates")

                # Gas
                gas_S.append(f['gas/Sigma'][:])
                gas_T.append(f['gas/T'][:])
                gas_cs.append(f['gas/cs'][:])
                gas_eta.append(f['gas/eta'][:])
                gas_nu.append(f['gas/nu'][:])

                # Snowlines (escalares por snapshot)
                for sp in ('H2O', 'CO2', 'CO'):
                    key = f'grid/rsnow_{sp}'
                    rsnow[sp].append(float(f[key][()]) if key in f else np.nan)

                # Dust
                dust_Sigma.append(f['dust/Sigma'][:])     # (Nr, 2)
                dust_vr.append(f['dust/v/rad'][:])         # (Nr, 5)
                dust_St.append(f['dust/St'][:])            # (Nr, 5)

                # ── Componentes químicos ───────────────────────────────────────────
                # Detectar especies volátiles disponibles en el primer archivo
                if r_grid is not None and not comp_sigma:
                    grid_keys = list(f['grid'].keys())
                    detected  = [k.replace('SigmaDust_', '') for k in grid_keys
                                 if k.startswith('SigmaDust_')]
                    vol_sps   = detected if detected else ['H2O', 'CO2', 'CO']
                    comp_sigma = {sp: [] for sp in vol_sps}
                    comp_sigma['silicates'] = []

                dust_snap    = f['dust/Sigma'][:].sum(-1)   # (Nr,) total polvo
                vol_total    = np.zeros_like(dust_snap)

                for sp in [k for k in comp_sigma if k != 'silicates']:
                    new_key = f'grid/SigmaDust_{sp}'
                    ice_key = f'grid/SigmaIce_{sp}'
                    old_key = f'components/{sp}/dust/Sigma'
                    if new_key in f:
                        arr = f[new_key][:]
                    elif ice_key in f:
                        arr = f[ice_key][:]
                    elif old_key in f:
                        arr = f[old_key][:]
                    else:
                        arr = np.zeros(len(r_grid))
                    comp_sigma[sp].append(arr)
                    vol_total += arr

                # Silicatos como remanente del polvo total
                comp_sigma['silicates'].append(np.maximum(0.0, dust_snap - vol_total))

        Nt = len(times_list)
        Nr = len(r_grid)
        print(f"[from_datadir] {Nt} snapshots válidos  |  Nr = {Nr}")

        # ── Construir instancia sin llamar __init__ ─────────────────────
        obj                 = cls.__new__(cls)
        obj.r               = r_grid
        obj.times           = np.array(times_list)
        obj.Nt              = Nt
        obj.M_star          = M_star * cls.M_SUN

        for sp in ('H2O', 'CO2', 'CO'):
            rsnow[sp] = np.array(rsnow[sp])
        obj.rsnow = rsnow

        obj.gas = {
            'Sigma': np.array(gas_S),    # (Nt, Nr)
            'T':     np.array(gas_T),
            'cs':    np.array(gas_cs),
            'eta':   np.array(gas_eta),
            'nu':    np.array(gas_nu),
        }
        obj.dust = {
            'Sigma': np.array(dust_Sigma),   # (Nt, Nr, 2)
            'vr':    np.array(dust_vr),       # (Nt, Nr, 5)
            'St':    np.array(dust_St),       # (Nt, Nr, 5)
        }

        # Componentes: verificar si se guardó algún Sigma
        if comp_sigma and all(len(v) > 0 for v in comp_sigma.values()):
            obj.comp = {sp: np.stack(comp_sigma[sp]) for sp in comp_sigma}
            obj._has_comp_sigma = True
            vol_sps = [s for s in obj.comp if s != 'silicates']
            print(f"[from_datadir] ✅ Componentes cargados: {vol_sps} + silicates")
        else:
            obj.comp = {}
            obj._has_comp_sigma = False
            print("[from_datadir] ⚠️  Componentes dust.Sigma NO guardados.")
            print("               → Composición estimada desde snowlines + abundancias iniciales.")

        # Derivados: Omega_K, H_gas, alpha
        obj.Omega_K = np.sqrt(obj.G_CGS * obj.M_star / obj.r**3)
        obj.H_gas   = obj.gas['cs'] / obj.Omega_K[np.newaxis, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            alpha_raw = obj.gas['nu'] / (obj.gas['cs'] * obj.H_gas)
        obj.alpha = np.clip(np.where(np.isfinite(alpha_raw), alpha_raw, 1e-3),
                            1e-5, 1e-1)
        return obj

    def _comp_from_snowlines(self, t_idx, r_emb):
        """
        Estima las fracciones de composición cuando component.dust.Sigma
        no está disponible en el HDF5.

        Lógica: si el embrión está FUERA del snowline de la especie X
        (r_emb > rsnow_X), la especie X está sólida y contribuye al pebble.
        Si está DENTRO (r_emb < rsnow_X), el hielo se sublimó → f_X = 0.
        Los silicatos siempre son refractarios (f_sil > 0 siempre).

        Las fracciones se normalizan a 1.
        """
        fracs = {}
        for sp, abun in self._ABUNDANCES.items():
            if sp == 'CO':      # CO no está en self.comp por defecto
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
    # Helpers internos
    # ══════════════════════════════════════════════════════════════════════

    def _interp(self, field_1d, r_emb):
        """Interpola en log(r) para suavidad. Devuelve float."""
        return float(np.interp(np.log(r_emb), np.log(self.r), field_1d))

    def _local(self, t, r_emb):
        """Devuelve dict con parámetros locales del disco en (t, r_emb)."""
        I = lambda arr: self._interp(arr[t], r_emb)
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, 1], r_emb)
        cs    = I(self.gas['cs'])
        eta   = I(self.gas['eta'])
        alpha = self._interp(self.alpha[t], r_emb)
        St    = self._interp(self.dust['St'][t, :, -1], r_emb)
        H_gas = self._interp(self.H_gas[t], r_emb)
        Omega = float(np.interp(np.log(r_emb), np.log(self.r), self.Omega_K))
        v_K   = Omega * r_emb
        # delta_v = η v_K  [L&J 2012, Eq. 16] — headwind pebble-embrión
        delta_v = eta * v_K + 1e-5 * cs   # floor numérico
        return dict(Sigma_peb=Sigma_peb, cs=cs, eta=eta, alpha=alpha,
                    St=St, H_gas=H_gas, Omega=Omega, v_K=v_K, delta_v=delta_v)

    # ══════════════════════════════════════════════════════════════════════
    # Física central
    # ══════════════════════════════════════════════════════════════════════

    def _pebble_flux(self, t, r_emb):
        """
        Ṁ_peb = 2π r Σ_peb |v_r|  [g/s]
        Referencia: Lambrechts & Johansen (2012), A&A 544, A32.
        Usa el bin de polvo grande (índice 1) y v del mayor Stokes (índice -1).
        """
        Sigma_peb = self._interp(self.dust['Sigma'][t, :, 1], r_emb)
        v_r       = self._interp(self.dust['vr'][t, :, -1],   r_emb)
        return 2 * np.pi * r_emb * Sigma_peb * abs(v_r)

    def _accretion_rate(self, M_core, r_emb, t):
        """
        Ṁ_core [g/s] — régimen de settling con transición suave 2D↔3D.
        
        Formation, Evolution, and Dynamics of Young Solar Systems
        Martin Pessah, Oliver Gressel
        Seccion 7.2.3 

        Física:
          R_H     = r (M/3M*)^(1/3)                     [radio de Hill]
          b_Bondi = √St · GM/Δv²    [Ormel & Klahr 2010, Eq. A9]
          b_Hill  = (St/0.1)^(1/3) R_H  [L&J 2012, Eq. A10]
          R_acc   = min(R_H, b_Bondi, b_Hill)
          H_peb   = H_gas √(α/(α+St))   [Youdin & Lithwick 2007, Eq.12]
          ρ_peb   = Σ_peb / (√(2π) H_peb)
          Ṁ_2D   = 2 R_acc Σ_peb Δv
          Ṁ_3D   = π R_acc² ρ_peb Δv
          Transición: Ṁ = erf(z)Ṁ_2D + (1-erf(z))Ṁ_3D, z=R_acc/(√2 H_peb)
          [Liu & Ormel 2018, A&A 615, A178, Appendix A.1]
        """
        if M_core <= 0:
            return 0.0
        p  = self._local(t, r_emb)
        if p['Sigma_peb'] < 1e-30:
            return 0.0

        G, M, M_s = self.G_CGS, M_core, self.M_star
        St, dv, Sigma = p['St'], p['delta_v'], p['Sigma_peb']

        # ── Radios característicos ────────────────────────────────────────
        R_H         = r_emb * (M / (3 * M_s))**(1/3)
        R_acc_Bondi = np.sqrt(St) * G * M / dv**2          # O&K 2010 Eq. A9
        R_acc_Hill  = (St / 0.1)**(1/3) * R_H              # L&J 2012 Eq. A10
        R_acc       = min(R_H, R_acc_Bondi, R_acc_Hill)

        # ── Escala vertical del disco de pebbles ──────────────────────────
        # Youdin & Lithwick (2007), Eq. 12
        H_peb = p['H_gas'] * np.sqrt(p['alpha'] / (p['alpha'] + St))
        H_peb = max(H_peb, 1e-10 * p['H_gas'])

        rho_peb = Sigma / (np.sqrt(2 * np.pi) * H_peb)

        # ── Tasas 2D y 3D ────────────────────────────────────────────────
        Mdot_2D = 2 * R_acc * Sigma * dv
        Mdot_3D = np.pi * R_acc**2 * rho_peb * dv

        # ── Transición suave: Liu & Ormel (2018), Appendix A.1 ───────────
        z   = R_acc / (np.sqrt(2) * H_peb)
        f2D = float(erf(z))
        return max(f2D * Mdot_2D + (1 - f2D) * Mdot_3D, 0.0)

    def _isolation_mass(self, r_emb, t):
        """
        Masa de aislamiento M_iso [g].

        Bitsch et al. (2018), A&A 612, A30, Eq. 5:
          M_iso = 25 M⊕ (h/0.05)³ · f_α · f_P
          f_α  = 0.34 [log₁₀(α/10⁻³) + 3] + 0.66
          f_P  = 1 − (d lnP/d lnr + 2.5) / 6
          d lnP/d lnr = −2η / h²   (de la definición de η)

        Fallback inferior: Lambrechts et al. (2014), A&A 572, A35:
          M_iso = 20 M⊕ (h/0.05)³
        """
        p     = self._local(t, r_emb)
        h     = p['H_gas'] / r_emb               # aspect ratio H/r
        alpha = p['alpha']
        eta   = p['eta']

        # Factor de viscosidad — Bitsch+2018 Eq. 5
        f_alpha    = 0.34 * (np.log10(np.clip(alpha, 1e-5, 1) / 1e-3) + 3) + 0.66

        # Factor de gradiente de presión
        d_lnP      = -2 * eta / (h**2 + 1e-30)    # d ln P / d ln r
        f_pressure = 1.0 - (d_lnP + 2.5) / 6.0

        M_Bitsch   = 25 * (h / 0.05)**3 * f_alpha * f_pressure * self.M_EARTH
        M_Lam14    = 20 * (h / 0.05)**3 * self.M_EARTH   # Lambrechts+2014 fallback
        return max(M_Bitsch, M_Lam14, 0.1 * self.M_EARTH)

    # ══════════════════════════════════════════════════════════════════════
    # API pública
    # ══════════════════════════════════════════════════════════════════════

    def run_growth(self, embryo_locations_AU, M0_g=1e24):
        """
        Simula el crecimiento simultáneo de todos los embriones sobre el disco
        pregrabado de tripodpy, con flujo de pebbles compartido correctamente.

        Lógica de flujo compartido (Ormel & Klahr 2010; Liu & Ormel 2018):
          - Los pebbles derivan RADIALMENTE HACIA ADENTRO.
          - Un embrión en r_out intercepta pebbles ANTES de que lleguen a r_in < r_out.
          - El flujo disponible para un embrión = Ṁ_disco(r) - Σ consumo de
            todos los embriones en r' > r en ese mismo paso de tiempo.
          - El loop anterior (embrión por embrión) usaba el mismo Ṁ_disco para
            TODOS, sobreestimando la masa acretada por los embriones interiores.

        Estructura:
          LOOP externo → snapshots i (disco estático en cada paso)
            LOOP interno → embriones ordenados de afuera hacia adentro
              - flux_consumed  se acumula con cada embrión procesado
              - Mdot_avail(r) = max(0, Mdot_disco(r) - flux_consumed)

        Parameters
        ----------
        embryo_locations_AU : list of float — radios [AU].
        M0_g : float — masa inicial [g]. Default 1e24 g ≈ 0.17 M_luna.

        Returns
        -------
        results : dict {r_AU: ndarray(N, 7)}
          Columnas: [t_s, M_core_g, M_H2O_g, M_CO2_g, M_sil_g,
                     r_snow_H2O_AU, M_iso_g]
        """
        # ── Ordenar embriones de afuera hacia adentro ─────────────────────────
        # Los pebbles vienen del disco externo. El embrión más exterior consume
        # primero; los interiores solo ven el remanente.
        locs_outer_to_inner = sorted(embryo_locations_AU, reverse=True)

        # ── Estado compartido de todos los embriones ──────────────────────────
        M_core   = {r: float(M0_g)              for r in locs_outer_to_inner}
        M_X      = {r: {sp: 0.0 for sp in self.comp} for r in locs_outer_to_inner}
        active   = {r: True                      for r in locs_outer_to_inner}
        histories= {r: []                        for r in locs_outer_to_inner}

        # ── Loop principal: snapshot por snapshot ─────────────────────────────
        for i in range(self.Nt):
            dt = float(self.times[i] - (self.times[i-1] if i > 0 else 0.0))
            if dt <= 0:
                continue

            # Snowline H2O en este tiempo (para el historial)
            r_snow_AU = float(self.rsnow['H2O'][i]) / self.AU \
                        if not np.isnan(self.rsnow.get('H2O', [np.nan])[i]) else np.nan

            # Consumo acumulado de pebbles [g/s] por embriones YA procesados
            # (siempre mayor radio → menor radio, afuera → adentro)
            flux_consumed = 0.0

            for r_au in locs_outer_to_inner:
                if not active[r_au]:
                    continue

                r_emb = r_au * self.AU

                # ── Masa de aislamiento ────────────────────────────────────────
                M_iso = self._isolation_mass(r_emb, i)
                if M_core[r_au] >= M_iso:
                    active[r_au] = False
                    continue

                # ── Flujo disponible: disco - lo que consumieron los exteriores ─
                Mdot_peb_disk  = self._pebble_flux(i, r_emb)
                Mdot_peb_avail = max(0.0, Mdot_peb_disk - flux_consumed)

                # ── Tasa de acreción ───────────────────────────────────────────
                Mdot_core_r = self._accretion_rate(M_core[r_au], r_emb, i)
                Mdot_core_r = min(Mdot_core_r, Mdot_peb_avail)  # limitado por flux real

                dM = Mdot_core_r * dt

                # ── Clamp: M_core no puede superar M_iso en un solo paso ───────
                # Sin este límite, un dt grande puede acumular masa muy por encima
                # de M_iso antes de que el check al inicio del loop lo detecte.
                # Esto hace la integración independiente de la resolución temporal.
                dM = min(dM, max(0.0, M_iso - M_core[r_au]))

                # Actualizar consumo acumulado para el siguiente embrión (interior)
                flux_consumed += Mdot_core_r

                # ── Composición química del dM acretado ───────────────────────
                if self._has_comp_sigma:
                    sig_ice_total = sum(
                        self._interp(self.comp[sp][i], r_emb) for sp in self.comp
                    )
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
                    self.times[i],
                    M_core[r_au],
                    M_X[r_au].get('H2O', 0),
                    M_X[r_au].get('CO2', 0),
                    M_X[r_au].get('silicates', 0),
                    r_snow_AU,
                    M_iso,
                ))

        # ── Reconstruir resultados en el orden original pedido ────────────────
        results = {
            r_au: (np.array(histories[r_au]) if histories[r_au] else np.empty((0, 7)))
            for r_au in embryo_locations_AU
        }
        return results


    def summary(self, results):
        """Imprime tabla resumen de composición final."""
        print(f"\n{'─'*70}")
        print(f"{'r [AU]':>8} {'M [M⊕]':>9} {'H₂O%':>7} {'CO₂%':>7} "
              f"{'Sil%':>7}  Tipo")
        print(f"{'─'*70}")
        for r_au, hist in results.items():
            if len(hist) == 0:
                print(f"{r_au:>8.2f}  — sin acreción")
                continue
            _, M, H2O, CO2, sil, _, _ = hist[-1]
            f_H2O = 100 * H2O / (M + 1e-30)
            f_CO2 = 100 * CO2 / (M + 1e-30)
            f_sil = 100 * sil / (M + 1e-30)
            tipo  = "🌊 Waterworld" if f_H2O > 10 else "🪨 Rocoso"
            print(f"{r_au:>8.2f}  {M/self.M_EARTH:>8.3f}  {f_H2O:>6.1f}  "
                  f"{f_CO2:>6.1f}  {f_sil:>6.1f}  {tipo}")
        print(f"{'─'*70}\n")
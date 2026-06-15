import numpy as np
import h5py
import dustpy.constants as c


class PebbleAccretionModule:

    def __init__(self, h5file, M_star=1.0):
        f = h5py.File(h5file, 'r')

        # =========================
        # GRID
        # =========================
        self.r = f['grid/r'][:]         # [cm]
        self.times = f['t'][:]          # [s]

        # Snowlines
        self.rsnow = {
            spec: f[f'grid/rsnow_{spec}'][:]
            for spec in ('H2O', 'CO2', 'CO')
        }

        # =========================
        # GAS
        # =========================
        self.gas = {
            key: f[f'gas/{key}'][:]
            for key in ('Sigma', 'T', 'cs', 'eta', 'nu')
        }

        # =========================
        # DUST
        # =========================
        self.dust = {
            'Sigma': f['dust/Sigma'][:],
            'vr': f['dust/v.rad'][:],
            'St': f['dust/St'][:]
        }

        # =========================
        # COMPOSICIÓN
        # =========================
        self.comp = {
            spec: f[f'components/{spec}/dust/Sigma'][:]
            for spec in ('H2O', 'CO2', 'silicates')
        }

        f.close()

        # =========================
        # ESTRELLA
        # =========================
        self.M_star = M_star * c.M_sun

        # =========================
        # DERIVADOS DEL DISCO
        # =========================
        self.Omega = np.sqrt(c.G * self.M_star / self.r**3)
        self.vK = self.r * self.Omega

    # =========================
    # INTERPOLADOR LOCAL
    # =========================
    def interp(self, arr, t_idx, r_emb):
        return np.interp(r_emb, self.r, arr[t_idx])

    # =========================
    # RADIOS FÍSICOS
    # =========================
    def hill_radius(self, M, r):
        return r * (M / (3 * self.M_star))**(1/3)

    def bondi_radius(self, M, delta_v):
        return c.G * M / (delta_v**2 + 1e-20)

    # =========================
    # ALTURA DE PEBBLES
    # =========================
    def pebble_scale_height(self, t_idx, r_emb, St):
        cs = self.interp(self.gas['cs'], t_idx, r_emb)
        Omega = np.interp(r_emb, self.r, self.Omega)
        nu = self.interp(self.gas['nu'], t_idx, r_emb)

        H_g = cs / Omega
        alpha = np.clip(nu / (cs * H_g), 1e-5, 1e-1)

        return H_g * np.sqrt(alpha / (alpha + St))

    # =========================
    # FLUJO DE PEBBLES
    # =========================
    def pebble_flux(self, t_idx, r_emb):

        Sigma_peb = np.interp(
            r_emb, self.r, self.dust['Sigma'][t_idx, :, 1]
        )

        v_r = np.interp(
            r_emb, self.r, self.dust['vr'][t_idx, :, -1]
        )

        return 2 * np.pi * r_emb * Sigma_peb * abs(v_r)

    # =========================
    # VELOCIDAD RELATIVA
    # =========================
    def relative_velocity(self, t_idx, r_emb, e, inc):

        eta = self.interp(self.gas['eta'], t_idx, r_emb)
        vK = np.interp(r_emb, self.r, self.vK)

        return np.sqrt(
            (eta * vK)**2 +
            (e * vK)**2 +
            (inc * vK)**2
        )

    # =========================
    # TASA DE ACRECIÓN
    # =========================
    def accretion_rate(self, M_core, r_emb, t_idx, e=0.0, inc=0.0):

        Sigma_peb = self.interp(self.dust['Sigma'][:, :, 1], t_idx, r_emb)
        St = self.interp(self.dust['St'][:, :, -1], t_idx, r_emb)

        delta_v = self.relative_velocity(t_idx, r_emb, e, inc)

        # Radios
        R_H = self.hill_radius(M_core, r_emb)
        R_B = self.bondi_radius(M_core, delta_v)

        R_acc = min(R_H, R_B * np.sqrt(St))

        # Altura pebbles
        H_peb = self.pebble_scale_height(t_idx, r_emb, St)

        # =========================
        # TRANSICIÓN SUAVE 2D-3D (Liu & Ormel-like)
        # =========================
        z = R_acc / (np.sqrt(2) * H_peb + 1e-20)
        f2D = np.tanh(z)  # más estable que erf

        # 2D
        Mdot_2D = 2 * R_acc * Sigma_peb * delta_v

        # 3D
        rho_peb = Sigma_peb / (np.sqrt(2*np.pi) * H_peb + 1e-20)
        Mdot_3D = np.pi * R_acc**2 * rho_peb * delta_v

        Mdot = f2D * Mdot_2D + (1 - f2D) * Mdot_3D

        return Mdot

    # =========================
    # MASA DE AISLAMIENTO (Bitsch+18)
    # =========================
    def pebble_isolation_mass(self, t_idx, r_emb):

        cs = self.interp(self.gas['cs'], t_idx, r_emb)
        Omega = np.interp(r_emb, self.r, self.Omega)

        H = cs / Omega
        h = H / r_emb

        return 20 * (h / 0.05)**3 * c.M_earth

    # =========================
    # CRECIMIENTO
    # =========================
    def run_growth(self, embryos_au, M0=1e24):

        results = {}

        for r_au in embryos_au:

            r_emb = r_au * c.au

            M_core = M0
            M_X = {spec: 0.0 for spec in self.comp}

            history = []

            for i, t in enumerate(self.times):

                Mdot_peb = self.pebble_flux(i, r_emb)
                Mdot_core = self.accretion_rate(M_core, r_emb, i)

                Mdot_core = min(Mdot_core, Mdot_peb)

                # Aislamiento
                M_iso = self.pebble_isolation_mass(i, r_emb)
                if M_core >= M_iso:
                    break

                dt = self.times[i] - (self.times[i-1] if i > 0 else 0)

                dM = Mdot_core * dt
                M_core += dM

                # =========================
                # COMPOSICIÓN LOCAL
                # =========================
                Sigma_tot = self.interp(
                    self.dust['Sigma'][:, :, 1], i, r_emb
                )

                if Sigma_tot < 1e-20:
                    continue

                for spec in self.comp:
                    frac = self.interp(self.comp[spec], i, r_emb) / (Sigma_tot + 1e-20)
                    M_X[spec] += frac * dM

                history.append([
                    t, M_core,
                    M_X['H2O'],
                    M_X['CO2'],
                    M_X['silicates']
                ])

            results[r_au] = np.array(history)

        return results
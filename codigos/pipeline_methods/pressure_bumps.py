import numpy as np
import dustpy.constants as c
from dustpylib.substructures.gaps import kanagawa2017, duffell2020

"""
Referencias
-----------
  Kanagawa et al. (2017) ApJ 838, 15
  Duffell (2020) ApJL 889, L11
  dustpylib docs: https://dustpylib.readthedocs.io/en/latest/
"""

class PressureBumpsMixin:

    # ════════════════════════════════════════════════════════════════════════
    # Helpers internos
    # ════════════════════════════════════════════════════════════════════════

    def _compute_gap_profile(self, gap_model, _M_p, _a, alpha0):

        q_now = _M_p / float(self.sim.star.M)

        h_p = float(
            np.interp(
                _a,
                self.sim.grid.r,
                self.sim.gas.Hp / self.sim.grid.r
            )
        )

        alp_p = float(
            np.interp(
                _a,
                self.sim.grid.r,
                alpha0
            )
        )

        f_gap = gap_model(
            self.sim.grid.r,
            _a,
            q_now,
            h_p,
            alp_p
        )

        f_safe = np.maximum(f_gap, 1e-3)

        return alpha0 / f_safe

    def _apply_static_alpha(self, alpha_static):

        self.sim.gas.alpha[...] = alpha_static
        self.sim.gas.alpha.updater = None

    def _imprint_gap(self, alpha_static, alpha0, label="gap"):

        ratio = alpha_static / alpha0

        print(f"  → Imprimiendo {label} en Σ_gas y Σ_dust...")

        self.sim.gas.Sigma[...]  /= ratio
        self.sim.dust.Sigma[...] /= ratio[:, None]

        self.sim.update()

    # ════════════════════════════════════════════════════════════════════════
    # Kanagawa et al. (2017)
    # ════════════════════════════════════════════════════════════════════════

    def setup_gap_kanagawa(
        self,
        M_planet:    float = None,
        a_planet_au: float = None,
        alpha_ref:   float = None,
        imprint:     bool  = False,
    ):

        _M_p = (
            M_planet
            if M_planet is not None
            else self.gap_M_planet
        )

        _a = (
            a_planet_au
            if a_planet_au is not None
            else self.gap_a_planet_au
        ) * c.au

        _a0_val = (
            alpha_ref
            if alpha_ref is not None
            else self.gap_alpha_ref
        )

        alpha0 = self.sim.gas.alpha.copy()

        M_star = float(self.sim.star.M)
        q      = _M_p / M_star

        print("Configurando gap planetario — Kanagawa et al. (2017)")
        print(
            f"  → M_planet  = {_M_p / c.M_earth:.2f} M⊕  "
            f"= {_M_p / c.M_jup:.4f} M_Jup"
        )
        print(f"  → a_planet  = {_a / c.au:.2f} AU")
        print(
            f"  → q = M_planet/M_star = {q:.3e}  "
            f"(M_star = {M_star / c.M_sun:.2f} M☉)"
        )
        print(
            f"  → alpha_ref = {_a0_val:.2e}  "
            f"[capturado en alpha0.copy()]"
        )

        alpha_static = self._compute_gap_profile(
            kanagawa2017,
            _M_p,
            _a,
            alpha0
        )

        self._apply_static_alpha(alpha_static)

        if imprint:

            self._imprint_gap(
                alpha_static,
                alpha0,
                label="gap Kanagawa"
            )

        alpha_max = float(self.sim.gas.alpha.max())

        print(
            f"  → α_max en gap = {alpha_max:.3e}  "
            f"|  α_ref = {_a0_val:.3e}"
        )

    # ════════════════════════════════════════════════════════════════════════
    # Duffell (2020)
    # ════════════════════════════════════════════════════════════════════════

    def setup_gap_duffell(
        self,
        M_planet:    float = None,
        a_planet_au: float = None,
        alpha_ref:   float = None,
        imprint:     bool  = False,
    ):

        _M_p = (
            M_planet
            if M_planet is not None
            else self.gap_M_planet
        )

        _a = (
            a_planet_au
            if a_planet_au is not None
            else self.gap_a_planet_au
        ) * c.au

        _a0_val = (
            alpha_ref
            if alpha_ref is not None
            else self.gap_alpha_ref
        )

        alpha0 = self.sim.gas.alpha.copy()

        M_star = float(self.sim.star.M)
        q      = _M_p / M_star

        print("Configurando gap planetario — Duffell (2020)")
        print(
            f"  → M_planet  = {_M_p / c.M_earth:.2f} M⊕  "
            f"= {_M_p / c.M_jup:.4f} M_Jup"
        )
        print(f"  → a_planet  = {_a / c.au:.2f} AU")
        print(
            f"  → q = M_planet/M_star = {q:.3e}  "
            f"(M_star = {M_star / c.M_sun:.2f} M☉)"
        )
        print(
            f"  → alpha_ref = {_a0_val:.2e}  "
            f"[capturado en alpha0.copy()]"
        )

        alpha_static = self._compute_gap_profile(
            duffell2020,
            _M_p,
            _a,
            alpha0
        )

        self._apply_static_alpha(alpha_static)

        if imprint:

            self._imprint_gap(
                alpha_static,
                alpha0,
                label="gap Duffell"
            )

        alpha_max = float(self.sim.gas.alpha.max())

        print(
            f"  → α_max en gap = {alpha_max:.3e}  "
            f"|  α_ref = {_a0_val:.3e}"
        )

    # ════════════════════════════════════════════════════════════════════════
    # Duffell múltiple
    # ════════════════════════════════════════════════════════════════════════

    def setup_gap_duffell_multi(
        self,
        planets: list,
        alpha_ref=None,
        imprint=False
    ):

        _a0 = (
            alpha_ref
            if alpha_ref is not None
            else self.gap_alpha_ref
        )

        print(
            f"Configurando {len(planets)} gaps planetarios "
            f"ESTÁTICOS — Duffell (multi)"
        )

        alpha0 = self.sim.gas.alpha.copy()

        f_safe_total = np.ones_like(self.sim.grid.r)

        for p in planets:

            Mp   = p["M_planet"]
            a_pl = p["a_planet_au"] * c.au

            alpha_partial = self._compute_gap_profile(
                duffell2020,
                Mp,
                a_pl,
                alpha0
            )

            f_safe_total *= alpha0 / alpha_partial

        alpha_static = alpha0 / f_safe_total

        self._apply_static_alpha(alpha_static)

        if imprint:

            self._imprint_gap(
                alpha_static,
                alpha0,
                label="gaps múltiples"
            )

    # ════════════════════════════════════════════════════════════════════════
    # Alpha sinusoidal
    # ════════════════════════════════════════════════════════════════════════

    def setup_alpha_sinusoidal(
        self,
        alpha_ref:   float = None,
        amplitude:   float = 5.0,
        n_bumps:     int   = 5,
        r_inner_au:  float = 5.0,
        r_outer_au:  float = 100.0,
        imprint:     bool  = False,
    ):

        _a0 = (
            alpha_ref
            if alpha_ref is not None
            else self.gap_alpha_ref
        )

        _r_in  = r_inner_au * c.au
        _r_out = r_outer_au * c.au
        _A     = amplitude
        _n     = int(n_bumps)

        print("Configurando α sinusoidal — gaps múltiples equidistantes")
        print(f"  → alpha_ref  = {_a0:.2e}")
        print(
            f"  → amplitud A = {_A:.1f}  "
            f"→  α_max = {_a0 * (1 + _A):.2e}"
        )
        print(
            f"  → n_gaps     = {_n}  "
            f"en [{r_inner_au:.1f}, {r_outer_au:.1f}] AU"
        )
        print(f"  → imprint    = {imprint}")

        r = self.sim.grid.r

        alpha_old = self.sim.gas.alpha.copy()
        alpha_static = alpha_old.copy()

        mask = (r >= _r_in) & (r <= _r_out)

        r_zone = r[mask]

        if r_zone.size > 0:

            log_in  = np.log(_r_in)
            log_out = np.log(_r_out)

            x = (
                (np.log(r_zone) - log_in)
                / (log_out - log_in)
            )

            perturb = (
                1.0
                + _A * np.sin(_n * np.pi * x) ** 2
            )

            alpha_static[mask] = alpha_old[mask] * perturb

        self._apply_static_alpha(alpha_static)

        if imprint:

            ratio = alpha_static / alpha_old
            ratio = np.where(ratio > 0, ratio, 1.0)

            _sf_gas  = getattr(
                self.sim.gas,
                'SigmaFloor',
                1e-100
            )

            _sf_dust = getattr(
                self.sim.dust,
                'SigmaFloor',
                1e-100
            )

            sigma_floor_gas = float(
                np.atleast_1d(
                    np.asarray(_sf_gas, dtype=float)
                ).min()
            )

            sigma_floor_dust = float(
                np.atleast_1d(
                    np.asarray(_sf_dust, dtype=float)
                ).min()
            )

            print(
                "  → Imprimiendo perfil sinusoidal "
                "estático en Σ_gas y Σ_dust..."
            )

            new_sigma_gas = (
                self.sim.gas.Sigma / ratio
            )

            new_sigma_dust = (
                self.sim.dust.Sigma / ratio[:, None]
            )

            self.sim.gas.Sigma[...] = np.maximum(
                new_sigma_gas,
                sigma_floor_gas
            )

            self.sim.dust.Sigma[...] = np.maximum(
                new_sigma_dust,
                sigma_floor_dust
            )

            self.sim.update()

        alpha_max = float(self.sim.gas.alpha.max())
        alpha_min = float(self.sim.gas.alpha.min())

        print(
            f"  → α ∈ [{alpha_min:.2e}, {alpha_max:.2e}]  "
            f"|  α_ref = {_a0:.2e}"
        )

    # ════════════════════════════════════════════════════════════════════════
    # Reset
    # ════════════════════════════════════════════════════════════════════════

    def reset_gap(self):

        _a0 = self.gap_alpha_ref

        print(
            f"Restaurando gas.alpha uniforme = {_a0:.2e} "
            f"(sin gap)..."
        )

        self.sim.gas.alpha[...] = np.full_like(
            self.sim.grid.r,
            _a0
        )

        self.sim.gas.alpha.updater = None

        self.sim.update()

        print("  → gap eliminado.")
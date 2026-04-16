"""
pressure_bumps.py — Gaps planetarios via modificación de gas.alpha
===================================================================

Mixin para WaterworldPipeline que implementa la formación de gaps planetarios
modificando el perfil de viscosidad del gas (gas.alpha updater).

Patrón físico (dustpylib docs oficiales)
-----------------------------------------
En estado estacionario: ν · Σ = cte  →  Σ ∝ 1/α

Para imponer un gap de profundidad f(r) = Σ_gap/Σ_0 se necesita:
    α(r) = α_0(r) / f(r)

El updater:
  1. Parte del perfil NO perturbado: alpha = alpha0.copy()
  2. Interpola h = H_gas/r  y  alp = alpha0  en la posición del planeta
     (escalares, NO arrays — así lo esperan kanagawa2017/duffell2020)
  3. Divide por el perfil de gap: alpha /= kanagawa2017(r, a, q, h, alp)

Opcionalmente se puede imprimir el gap en Σ_gas y Σ_dust desde t=0
(en lugar de esperar al tiempo viscoso).

Métodos disponibles
-------------------
  setup_gap_kanagawa(M_planet, a_planet_au, alpha_ref, imprint)
  setup_gap_duffell(M_planet, a_planet_au, alpha_ref, imprint)
  reset_gap()

M_planet en CGS. Ejemplos de uso:
    pipeline.setup_gap_kanagawa(M_planet=1.0 * c.M_jup, a_planet_au=5.2)
    pipeline.setup_gap_duffell(M_planet=30. * c.M_earth, a_planet_au=5.)
    pipeline.setup_gap_kanagawa(M_planet=2. * c.M_earth, a_planet_au=3.)

Referencias
-----------
  Kanagawa et al. (2017) ApJ 838, 15
  Duffell (2020) ApJL 889, L11
  dustpylib docs: https://dustpylib.readthedocs.io/en/latest/
"""

import numpy as np
import dustpy.constants as c
from scipy.interpolate import interp1d
from dustpylib.substructures.gaps import kanagawa2017, duffell2020


class PressureBumpsMixin:
    """
    Mixin de gaps planetarios por modificación de α(r).
    Asume que self.sim está inicializado (initialize_simulation llamado).

    Parámetros de gap configurables en el __init__ del pipeline:
        gap_alpha_ref   : float  — α de referencia [adim], default 1e-3
        gap_a_planet_au : float  — semieje del planeta [AU], default 5.0
        gap_M_planet    : float  — masa del planeta [g], default 30*c.M_earth
    """

    # ══════════════════════════════════════════════════════════════════════════
    # Kanagawa et al. (2017)
    # ══════════════════════════════════════════════════════════════════════════

    def setup_gap_kanagawa(
        self,
        M_planet:    float = None,
        a_planet_au: float = None,
        alpha_ref:   float = None,
        imprint:     bool  = False,
    ):
        """
        Instala un gap planetario usando el modelo de Kanagawa et al. (2017).

        Sigue exactamente el patrón de la documentación oficial de dustpylib:
          - alpha0 se copia ANTES de instalar el updater
          - h y alp se INTERPOLAN en la posición del planeta (escalares)
          - El updater divide alpha0 element-wise por el perfil del gap

        Parameters
        ----------
        M_planet : float, optional
            Masa del planeta [g — CGS]. Ejemplos:
              1.0 * c.M_jup    → Júpiter
              0.3 * c.M_jup    → Saturno
              30.  * c.M_earth → ~0.1 M_Jup (pressure bump parcial)
            Default: self.gap_M_planet.
        a_planet_au : float, optional
            Semieje del planeta [AU]. Default: self.gap_a_planet_au.
        alpha_ref : float, optional
            α de referencia fuera del gap [adim]. Default: self.gap_alpha_ref.
            ⚠ Se captura en alpha0 ANTES del updater — no usa sim.gas.alpha
            dentro del closure para evitar circularidad.
        imprint : bool, optional
            Si True, aplica el gap a Σ_gas y Σ_dust desde t=0 dividiendo
            por alpha/alpha0 (en lugar de esperar al tiempo viscoso).
            Default: False.
        """
        _M_p = M_planet    if M_planet    is not None else self.gap_M_planet
        _a   = (a_planet_au if a_planet_au is not None
                else self.gap_a_planet_au) * c.au
        _a0_val = alpha_ref if alpha_ref is not None else self.gap_alpha_ref

        # ── Copiar alpha0 ANTES de instalar el updater (patrón doc oficial) ──
        # Puede ser array espacialmente variante; para alpha constante es
        # equivalente a un escalar, pero capturamos el array completo.
        alpha0 = self.sim.gas.alpha.copy()

        # Calcular q con la masa estelar actual
        M_star = float(self.sim.star.M)
        q      = _M_p / M_star

        # Diagnóstico
        print("Configurando gap planetario — Kanagawa et al. (2017)")
        print(f"  → M_planet  = {_M_p / c.M_earth:.2f} M⊕  "
              f"= {_M_p / c.M_jup:.4f} M_Jup")
        print(f"  → a_planet  = {_a / c.au:.2f} AU")
        print(f"  → q = M_planet/M_star = {q:.3e}  "
              f"(M_star = {M_star / c.M_sun:.2f} M☉)")
        print(f"  → alpha_ref = {_a0_val:.2e}  [capturado en alpha0.copy()]")

        # Capturar en el closure (alpha0 es el array unperturbado)
        _a_cl     = _a
        _alpha0   = alpha0    # array completo, constante en el closure

        def _gap_updater_kanagawa(sim):
            """
            Updater de sim.gas.alpha para el gap de Kanagawa 2017.

            Sigue el patrón oficial:
              1. Partir de alpha_unperturbado = alpha0.copy()
              2. Interpolar h y alp en la posición del planeta (escalares)
              3. alpha /= kanagawa2017(r, a, q, h_scalar, alp_scalar)
            """
            alpha_out = _alpha0.copy()

            q_now = _M_p / float(sim.star.M)   # dinámico con la estrella

            # Interpolar h = H_gas/r en la posición del planeta (escalar)
            f_h   = interp1d(sim.grid.r, sim.gas.Hp / sim.grid.r,
                             bounds_error=False, fill_value="extrapolate")
            h_p   = float(f_h(_a_cl))

            # Interpolar alpha0 en la posición del planeta (escalar)
            f_alp = interp1d(sim.grid.r, _alpha0,
                             bounds_error=False, fill_value="extrapolate")
            alp_p = float(f_alp(_a_cl))

            # Dividir por el perfil de gap (array de shape Nr)
            f_gap  = kanagawa2017(sim.grid.r, _a_cl, q_now, h_p, alp_p)
            f_safe = np.maximum(f_gap, 1e-10)
            alpha_out /= f_safe
            return alpha_out

        self.sim.gas.alpha.updater.updater = _gap_updater_kanagawa
        self.sim.gas.alpha.update()
        self.sim.update()

        # Imprimir gap en Σ desde t=0 (opcional)
        if imprint:
            ratio = self.sim.gas.alpha / alpha0
            print("  → Imprimiendo gap en Σ_gas y Σ_dust (imprint=True)...")
            self.sim.gas.Sigma[...]  /= ratio
            self.sim.dust.Sigma[...] /= ratio[:, None]
            self.sim.update()

        alpha_max = float(self.sim.gas.alpha.max())
        print(f"  → α_max en gap = {alpha_max:.3e}  |  α_ref = {_a0_val:.3e}")

    # ══════════════════════════════════════════════════════════════════════════
    # Duffell (2020)
    # ══════════════════════════════════════════════════════════════════════════

    def setup_gap_duffell(
        self,
        M_planet:    float = None,
        a_planet_au: float = None,
        alpha_ref:   float = None,
        imprint:     bool  = False,
    ):
        """
        Instala un gap planetario usando el modelo de Duffell (2020).

        Mismo patrón que setup_gap_kanagawa() pero usando duffell2020().
        Alternativa con una parametrización empírica diferente del perfil.

        Parameters
        ----------
        M_planet : float, optional
            Masa del planeta [g — CGS]. Default: self.gap_M_planet.
        a_planet_au : float, optional
            Semieje del planeta [AU]. Default: self.gap_a_planet_au.
        alpha_ref : float, optional
            α de referencia [adim]. Default: self.gap_alpha_ref.
        imprint : bool, optional
            Si True, aplica el gap a Σ desde t=0. Default: False.
        """
        _M_p    = M_planet    if M_planet    is not None else self.gap_M_planet
        _a      = (a_planet_au if a_planet_au is not None
                   else self.gap_a_planet_au) * c.au
        _a0_val = alpha_ref   if alpha_ref   is not None else self.gap_alpha_ref

        # Copiar alpha0 ANTES de instalar el updater
        alpha0 = self.sim.gas.alpha.copy()

        M_star = float(self.sim.star.M)
        q      = _M_p / M_star

        print("Configurando gap planetario — Duffell (2020)")
        print(f"  → M_planet  = {_M_p / c.M_earth:.2f} M⊕  "
              f"= {_M_p / c.M_jup:.4f} M_Jup")
        print(f"  → a_planet  = {_a / c.au:.2f} AU")
        print(f"  → q = M_planet/M_star = {q:.3e}  "
              f"(M_star = {M_star / c.M_sun:.2f} M☉)")
        print(f"  → alpha_ref = {_a0_val:.2e}  [capturado en alpha0.copy()]")

        _a_cl   = _a
        _alpha0 = alpha0

        def _gap_updater_duffell(sim):
            """
            Updater de sim.gas.alpha para el gap de Duffell 2020.

            Sigue el patrón oficial:
              1. Partir de alpha0.copy()
              2. Interpolar h y alp en la posición del planeta (escalares)
              3. alpha /= duffell2020(r, a, q, h_scalar, alp_scalar)
            """
            alpha_out = _alpha0.copy()

            q_now = _M_p / float(sim.star.M)

            f_h   = interp1d(sim.grid.r, sim.gas.Hp / sim.grid.r,
                             bounds_error=False, fill_value="extrapolate")
            h_p   = float(f_h(_a_cl))

            f_alp = interp1d(sim.grid.r, _alpha0,
                             bounds_error=False, fill_value="extrapolate")
            alp_p = float(f_alp(_a_cl))

            f_gap  = duffell2020(sim.grid.r, _a_cl, q_now, h_p, alp_p)
            f_safe = np.maximum(f_gap, 1e-10)
            alpha_out /= f_safe
            return alpha_out

        self.sim.gas.alpha.updater.updater = _gap_updater_duffell
        self.sim.gas.alpha.update()
        self.sim.update()

        if imprint:
            ratio = self.sim.gas.alpha / alpha0
            print("  → Imprimiendo gap en Σ_gas y Σ_dust (imprint=True)...")
            self.sim.gas.Sigma[...]  /= ratio
            self.sim.dust.Sigma[...] /= ratio[:, None]
            self.sim.update()

        alpha_max = float(self.sim.gas.alpha.max())
        print(f"  → α_max en gap = {alpha_max:.3e}  |  α_ref = {_a0_val:.3e}")

    # ══════════════════════════════════════════════════════════════════════════
    # Utilidades
    # ══════════════════════════════════════════════════════════════════════════

    def reset_gap(self):
        """
        Elimina el gap restaurando gas.alpha al perfil uniforme gap_alpha_ref.

        Útil para comparar runs con y sin gap, o cambiar de modelo sin
        reinicializar toda la simulación.
        """
        _a0 = self.gap_alpha_ref
        print(f"Restaurando gas.alpha uniforme = {_a0:.2e} (sin gap)...")
        self.sim.gas.alpha.updater.updater = lambda sim: np.full_like(sim.grid.r, _a0)
        self.sim.gas.alpha.update()
        self.sim.update()
        print("  → gap eliminado.")

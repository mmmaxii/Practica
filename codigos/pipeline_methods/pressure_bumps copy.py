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

        # Calcular gap de Kanagawa en t=0 de forma ESTÁTICA
        q_now = _M_p / float(self.sim.star.M)
        h_p = float(np.interp(_a_cl, self.sim.grid.r, self.sim.gas.Hp / self.sim.grid.r))
        alp_p = float(np.interp(_a_cl, self.sim.grid.r, alpha0))

        f_gap = kanagawa2017(self.sim.grid.r, _a_cl, q_now, h_p, alp_p)
        
        # Limitar la profundidad a max 99.9% para evitar 'Factor exactly singular'
        f_safe = np.maximum(f_gap, 1e-3)
        alpha_static = alpha0 / f_safe

        # Fijar vector y destruir updater continuo
        self.sim.gas.alpha[...] = alpha_static
        self.sim.gas.alpha.updater = None
        
        # IMPRINT FORZADO: tallar el disco en t=0 para evitar flujos viscosos violentos
        ratio = alpha_static / alpha0
        print("  → Imprimiendo gap en Σ_gas y Σ_dust de inmediato (Static Imprint)...")
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

        # Calcular gap de Duffell en t=0 de forma ESTÁTICA
        q_now = _M_p / float(self.sim.star.M)
        h_p = float(np.interp(_a_cl, self.sim.grid.r, self.sim.gas.Hp / self.sim.grid.r))
        alp_p = float(np.interp(_a_cl, self.sim.grid.r, alpha0))

        f_gap = duffell2020(self.sim.grid.r, _a_cl, q_now, h_p, alp_p)
        
        # Limitar la profundidad a max 99.9% para evitar 'Factor exactly singular'
        f_safe = np.maximum(f_gap, 1e-3)
        alpha_static = alpha0 / f_safe

        # Fijar vector y destruir updater continuo
        self.sim.gas.alpha[...] = alpha_static
        self.sim.gas.alpha.updater = None

        # IMPRINT FORZADO: tallar el gap enseguida para no congelar dt
        ratio = alpha_static / alpha0
        print("  → Imprimiendo gap de Duffell en Σ_gas y Σ_dust de inmediato (Static)...")
        self.sim.gas.Sigma[...]  /= ratio
        self.sim.dust.Sigma[...] /= ratio[:, None]
        self.sim.update()

        alpha_max = float(self.sim.gas.alpha.max())
        print(f"  → α_max en gap = {alpha_max:.3e}  |  α_ref = {_a0_val:.3e}")

    def setup_gap_duffell_multi(self, planets: list, alpha_ref=None, imprint=False):
        _a0 = alpha_ref if alpha_ref is not None else self.gap_alpha_ref
        _pl_list = [(p["M_planet"], p["a_planet_au"] * c.au) for p in planets]
        
        # Calcular gaps de Duffell para el sistema enjambre/múltiple
        print(f"Configurando {len(_pl_list)} gaps planetarios ESTÁTICOS — Duffell (multi)")

        f_safe_total = np.ones_like(self.sim.grid.r)
        alpha0 = self.sim.gas.alpha.copy()

        for Mp, a_pl in _pl_list:
            q      = Mp / float(self.sim.star.M)
            h_p    = float(np.interp(a_pl, self.sim.grid.r, self.sim.gas.Hp / self.sim.grid.r))
            alp_p  = float(np.interp(a_pl, self.sim.grid.r, alpha0))
            f_gap  = duffell2020(self.sim.grid.r, a_pl, q, h_p, alp_p)
            f_safe_total *= np.maximum(f_gap, 1e-3)
        
        alpha_static = alpha0 / f_safe_total
        self.sim.gas.alpha[...] = alpha_static
        self.sim.gas.alpha.updater = None

        if imprint:
            ratio = alpha_static / alpha0
            print("  → Imprimiendo gaps múltiples en Σ_gas y Σ_dust...")
            self.sim.gas.Sigma[...] /= ratio
            self.sim.dust.Sigma[...] /= ratio[:, None]
            self.sim.update()

    # ══════════════════════════════════════════════════════════════════════════
    # Alpha sinusoidal — múltiples gaps equidistantes
    # ══════════════════════════════════════════════════════════════════════════

    def setup_alpha_sinusoidal(
        self,
        alpha_ref:   float = None,
        amplitude:   float = 5.0,
        n_bumps:     int   = 5,
        r_inner_au:  float = 5.0,
        r_outer_au:  float = 100.0,
        imprint:     bool  = False,
    ):
        """
        Genera múltiples gaps equidistantes mediante una perturbación sinusoidal
        sobre el perfil de viscosidad α(r).

        El perfil resultante es:
            α(r) = α_ref · (1 + A · sin²(n_bumps · π · x(r)))

        donde x(r) ∈ [0, 1] es la coordenada radial normalizada en escala log
        entre r_inner y r_outer.

        Física:
          - Los máximos de α (α_ref·(1+A)) corresponden a gaps (Σ ∝ 1/α).
          - Los mínimos de α (α_ref) corresponden a anillos / trampas de presión.
          - A=0   → sin estructura (α uniforme)
          - A=1   → amplitud suave  (~2× el nivel base)
          - A=5   → amplitud media  (~6× el nivel base)
          - A=10  → amplitud fuerte (~11× el nivel base)

        Parameters
        ----------
        alpha_ref : float, optional
            α de referencia (nivel base de viscosidad) [adim].
            Default: self.gap_alpha_ref.
            ⚠ No usar sim.gas.alpha — se captura antes del updater.
        amplitude : float, optional
            Amplitud A de la perturbación sinusoidal [adim].
            Ejemplos: 1.0 (suave), 5.0 (media), 10.0 (fuerte).
            Default: 5.0.
        n_bumps : int, optional
            Número de gaps (semiperíodos de la sinusoidal en [r_inner, r_outer]).
            Default: 5.
        r_inner_au : float, optional
            Radio interior de la región perturbada [AU]. Default: 5 AU.
        r_outer_au : float, optional
            Radio exterior de la región perturbada [AU]. Default: 100 AU.
        imprint : bool, optional
            Si True, aplica el perfil α a Σ_gas y Σ_dust desde t=0.
            Default: False.

        Notes
        -----
        Fuera de [r_inner, r_outer] el α permanece en alpha_ref (sin estructura).
        La sinusoidal usa sin² para garantizar α ≥ α_ref siempre (sin valores
        negativos ni inferiores al baseline).
        """
        _a0        = alpha_ref  if alpha_ref  is not None else self.gap_alpha_ref
        _r_in      = r_inner_au * c.au
        _r_out     = r_outer_au * c.au
        _A         = amplitude
        _n         = int(n_bumps)

        # Capturar alpha0 array ANTES del updater (patrón dustpylib)
        alpha0 = self.sim.gas.alpha.copy()

        print("Configurando α sinusoidal — gaps múltiples equidistantes")
        print(f"  → alpha_ref  = {_a0:.2e}")
        print(f"  → amplitud A = {_A:.1f}  →  α_max = {_a0 * (1 + _A):.2e}")
        print(f"  → n_gaps     = {_n}  en [{r_inner_au:.1f}, {r_outer_au:.1f}] AU")
        print(f"  → imprint    = {imprint}")

        # Cerrar variables en el closure
        _a0_cl  = _a0
        _A_cl   = _A
        _n_cl   = _n
        _rin_cl = _r_in
        _rot_cl = _r_out

        # Calcular Sinusoidales
        r   = self.sim.grid.r
        alpha_static = np.full_like(r, _a0)
        mask   = (r >= _r_in) & (r <= _r_out)
        r_zone = r[mask]

        if r_zone.size > 0:
            log_in  = np.log(_r_in)
            log_out = np.log(_r_out)
            x       = (np.log(r_zone) - log_in) / (log_out - log_in)
            perturb = 1.0 + _A * np.sin(_n * np.pi * x) ** 2
            alpha_static[mask] = _a0 * perturb

        self.sim.gas.alpha[...] = alpha_static
        self.sim.gas.alpha.updater = None

        # Imprint forzado — con SigmaFloor para prevenir celdas nulas
        # que harían singular la matriz LU del solver implícito (RuntimeError).
        alpha_base = np.full_like(self.sim.grid.r, _a0)
        ratio = alpha_static / alpha_base
        ratio = np.where(ratio > 0, ratio, 1.0)

        # Leer SigmaFloor de la simulación (tripodpy lo define en gas y dust).
        # SigmaFloor puede ser un Field-array (Nr,) o un escalar según la versión
        # de tripodpy, por eso se extrae con np.atleast_1d().min() en vez de float().
        _sf_gas  = getattr(self.sim.gas,  'SigmaFloor', 1e-100)
        _sf_dust = getattr(self.sim.dust, 'SigmaFloor', 1e-100)
        sigma_floor_gas  = float(np.atleast_1d(np.asarray(_sf_gas,  dtype=float)).min())
        sigma_floor_dust = float(np.atleast_1d(np.asarray(_sf_dust, dtype=float)).min())

        print("  → Imprimiendo perfil sinusoidal estático en Σ_gas y Σ_dust...")
        new_sigma_gas  = self.sim.gas.Sigma  / ratio
        new_sigma_dust = self.sim.dust.Sigma / ratio[:, None]

        # Aplicar piso: nunca dejar celdas por debajo de SigmaFloor
        self.sim.gas.Sigma[...]  = np.maximum(new_sigma_gas,  sigma_floor_gas)
        self.sim.dust.Sigma[...] = np.maximum(new_sigma_dust, sigma_floor_dust)
        self.sim.update()

        alpha_max = float(self.sim.gas.alpha.max())
        alpha_min = float(self.sim.gas.alpha.min())
        print(f"  → α ∈ [{alpha_min:.2e}, {alpha_max:.2e}]  |  α_ref = {_a0:.2e}")
        print(f"  → Σ_gas  floor = {sigma_floor_gas:.2e}  |  Σ_dust floor = {sigma_floor_dust:.2e}")

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


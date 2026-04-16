"""
disk_setup.py — Configuración de grilla, estrella e inicialización de tripodpy
===============================================================================

Mixin para WaterworldPipeline que agrupa:
  - setup_grid()             → grilla radial
  - setup_star()             → parámetros estelares (M, R, T)
  - initialize_simulation()  → llama a sim.initialize()
  - setup_star_evolution()   → contracción pre-MS (Hayashi track)
"""

import numpy as np
import dustpy.constants as c


class DiskSetupMixin:
    """
    Mixin de configuración inicial del disco.
    Asume que self.sim es una instancia de tripodpy.Simulation
    ya creada en WaterworldPipeline.__init__().
    """

    # ══════════════════════════════════════════════════════════════════════════
    # Grilla radial
    # ══════════════════════════════════════════════════════════════════════════

    def setup_grid(self, rmin=1.0 * c.au, rmax=1000.0 * c.au, Nr=300):
        """
        Configura la grilla radial logarítmica.

        Parameters
        ----------
        rmin : float
            Radio interior [cm]. Default: 1 AU.
        rmax : float
            Radio exterior [cm]. Default: 1000 AU.
        Nr : int
            Número de celdas radiales. Default: 300.

        Notes
        -----
        Posteriormente se puede agregar refinamiento de grilla (dustpylib)
        si los gradientes en las snowlines se vuelven numéricamente difíciles.
        """
        print(f"Configurando grilla: {Nr} celdas  |  "
              f"{rmin/c.au:.1f} – {rmax/c.au:.1f} AU")
        self.sim.ini.grid.rmin = rmin
        self.sim.ini.grid.rmax = rmax
        self.sim.ini.grid.Nr   = Nr
        # Radio característico del perfil de Σ_gas (Lynden-Bell & Pringle 1974)
        self.sim.ini.gas.SigmaRc = 30.0 * c.au

    # ══════════════════════════════════════════════════════════════════════════
    # Estrella central
    # ══════════════════════════════════════════════════════════════════════════

    def setup_star(self, M_star_Msun=None, R_star_Rsun=None, T_star_K=None):
        """
        Inicializa los parámetros de la estrella central.

        Acepta valores en unidades solares / Kelvin y los convierte a CGS
        antes de asignarlos a tripodpy (que trabaja internamente en CGS).

        Parameters
        ----------
        M_star_Msun : float, optional
            Masa estelar [M_sun]. Default: self.M_star_Msun (1.0).
        R_star_Rsun : float, optional
            Radio estelar [R_sun]. Default: self.R_star_Rsun (2.0, radio joven
            inflado; Baraffe et al. 2015).
        T_star_K : float, optional
            Temperatura efectiva [K]. Default: self.T_star_K (5778 K).
        """
        if M_star_Msun is None: M_star_Msun = self.M_star_Msun
        if R_star_Rsun is None: R_star_Rsun = self.R_star_Rsun
        if T_star_K    is None: T_star_K    = self.T_star_K

        # Guardar en unidades solares para referencia posterior
        self.M_star_Msun = M_star_Msun
        self.R_star_Rsun = R_star_Rsun
        self.T_star_K    = T_star_K

        # Convertir a CGS y asignar a tripodpy
        self.sim.ini.star.M = M_star_Msun * c.M_sun
        self.sim.ini.star.R = R_star_Rsun * c.R_sun
        self.sim.ini.star.T = T_star_K

        print(f"Estrella configurada:  "
              f"M = {M_star_Msun:.2f} M☉  |  "
              f"R = {R_star_Rsun:.2f} R☉  |  "
              f"T = {T_star_K:.0f} K")

    # ══════════════════════════════════════════════════════════════════════════
    # Inicialización de tripodpy
    # ══════════════════════════════════════════════════════════════════════════

    def initialize_simulation(self):
        """
        Llama a sim.initialize() de tripodpy.

        DEBE ejecutarse después de setup_grid() y setup_star(), y ANTES de
        agregar componentes ad-hoc (add_volatile_components, etc.).
        """
        print("Inicializando la simulación base (tripodpy)...")
        self.sim.initialize()

    # ══════════════════════════════════════════════════════════════════════════
    # Evolución estelar (contracción pre-MS)
    # ══════════════════════════════════════════════════════════════════════════

    def setup_star_evolution(self, R_fin_Rsun=1.5, t_contract_yr=1.0e7):
        """
        Configura la contracción estelar pre-secuencia principal (Hayashi track).

        El radio estelar decrece linealmente desde R_ini (fijado en setup_star)
        hasta R_fin durante t_contract años. Pasado ese tiempo, R se fija en
        R_fin (llegada a la ZAMS).

        Parameters
        ----------
        R_fin_Rsun : float
            Radio final (ZAMS) [R_sun]. Default: 1.5.
        t_contract_yr : float
            Tiempo de Kelvin-Helmholtz [yr]. Default: 1e7 (Hayashi 1961).

        Notes
        -----
        Las constantes del closure se capturan UNA sola vez al llamar este
        método, no en cada timestep, evitando bugs si self.R_star_Rsun cambia
        posteriormente.
        """
        print("Configurando evolución estelar (contracción del radio)...")

        R_ini_cgs      = float(self.sim.star.R)          # [cm] — ya en CGS
        R_fin_cgs      = R_fin_Rsun    * c.R_sun          # [cm]
        t_contract_cgs = t_contract_yr * c.year           # [s]

        if R_ini_cgs <= R_fin_cgs:
            print(f"  [!] R_ini ({self.R_star_Rsun:.2f} R☉) ≤ R_fin ({R_fin_Rsun:.2f} R☉): "
                  f"no hay contracción posible. sim.star.R permanece estático.")
            return

        slope = (R_fin_cgs - R_ini_cgs) / t_contract_cgs   # < 0 (contracción)

        print(f"  → R_ini = {self.R_star_Rsun:.2f} R☉  →  R_fin = {R_fin_Rsun:.2f} R☉  "
              f"en {t_contract_yr:.1e} yr")
        print(f"  → dR/dt = {slope / c.R_sun * c.year:.3e} R☉/yr")

        def R_star_evolving(sim):
            R = R_ini_cgs + slope * sim.t
            return float(np.maximum(R, R_fin_cgs))   # clamp al radio ZAMS

        self.sim.star.R.updater.updater = R_star_evolving
        self.sim.star.R.update()

        r_now = float(self.sim.star.R) / c.R_sun
        print(f"  → R(t=0) = {r_now:.4f} R☉  (esperado: {self.R_star_Rsun:.2f} R☉)")

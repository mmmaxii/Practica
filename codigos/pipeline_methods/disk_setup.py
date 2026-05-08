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
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from grid_extension import refinegrid_dustpy



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
        print(f"Configurando grilla original (log): {Nr} celdas  |  "
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
    # Refinamiento de grilla
    # ══════════════════════════════════════════════════════════════════════════

    def setup_refined_grid(
        self,
        gap_positions_au: list = None,
        snowline_au: float = 2.0,
        num_gap: int = 3,
        num_snow: int = 2,
    ):
        """
        Refina la grilla radial antes de initialize_simulation().

        DEBE llamarse después de setup_grid()/setup_star() pero ANTES de
        initialize_simulation(). Sigue el patrón oficial de dustpy docs:

            sim.grid.ri = ri_refinada   # ndarray plano, ANTES de initialize()
            sim.initialize()            # makegrids() convierte ri a Field constante

        No modificar sim.ini.grid.Nr — makegrids() lo calcula automáticamente
        como len(ri) - 1.

        Parameters
        ----------
        gap_positions_au : list of float, optional
            Posiciones de los gaps [AU]. Se refina con num=num_gap en cada una.
        snowline_au : float
            Posición de la snowline H2O [AU] para refinar. Default: 2.0 AU.
        num_gap : int
            Niveles de refinamiento en cada gap. Default: 3.
        num_snow : int
            Niveles de refinamiento en la snowline. Default: 2.
        """
        from grid_extension import refinegrid_tripodpy   # función del tutorial dustpy

        rmin = float(self.sim.ini.grid.rmin)
        rmax = float(self.sim.ini.grid.rmax)
        Nr   = int(self.sim.ini.grid.Nr)

        # Construir ri log-uniforme (igual que haría initialize() por defecto)
        ri = np.geomspace(rmin, rmax, Nr + 1)

        # ── Refinamiento en snowline H2O ──────────────────────────────────
        ri = refinegrid_tripodpy(ri, snowline_au * c.au, num=num_snow)
        print(f"  → Grid refinado en snowline H2O ({snowline_au:.1f} AU)  num={num_snow}")

        # ── Refinamiento en cada posición de gap ──────────────────────────
        if gap_positions_au:
            for r_au in sorted(gap_positions_au):
                r_cm = r_au * c.au
                if rmin < r_cm < rmax:
                    ri = refinegrid_tripodpy(ri, r_cm, num=num_gap)
                    print(f"  → Grid refinado en gap ({r_au:.1f} AU)  num={num_gap}")
                else:
                    print(f"  [!] Gap en {r_au:.1f} AU fuera del dominio — omitido")

        Nr_new = len(ri) - 1

        # ── Asignar como ndarray plano ANTES de initialize() ─────────────
        # dustpy docs: "It is sufficient to assign a numpy.ndarray to
        # Simulation.grid.ri and not a simframe.Field."
        # makegrids() (dentro de initialize()) convierte ri a Field constante
        # y calcula Nr, r, A automáticamente desde len(ri)-1.
        self.sim.grid.ri = ri   # ndarray, NO tocar ini.grid.Nr

        print(f"  → Nr: {Nr} → {Nr_new} celdas  (+{Nr_new - Nr})")





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

    # ══════════════════════════════════════════════════════════════════════════
    # Reinicio desde dump (restart)
    # ══════════════════════════════════════════════════════════════════════════

    @staticmethod
    def restart_from_dump(
        datadir: str,
        t_end_years: float,
        num_extra_snapshots: int = 30,
        dump_filename: str = "frame.dmp",
    ):
        """
        Reinicia una simulación desde el dump file generado por tripodpy.

        El dump file contiene el estado completo de la simulación, incluyendo
        todos los updaters y funciones personalizadas. NO es necesario
        re-ejecutar ningún paso del pipeline.

        Patrón oficial (tripodpy docs):
            sim = readdump("frame.dmp")
            sim.t.snapshots = np.concatenate([sim.t.snapshots, nuevos_snaps])
            sim.run()

        Parameters
        ----------
        datadir : str
            Directorio de datos que contiene el dump file y los HDF5.
        t_end_years : float
            Nuevo tiempo final [yr]. Debe ser mayor al tiempo actual del dump.
        num_extra_snapshots : int
            Número de snapshots adicionales entre t_actual y t_end_years.
            Default: 30.
        dump_filename : str
            Nombre del archivo dump. Default: "frame.dmp".

        Returns
        -------
        sim : tripodpy.Simulation
            El objeto de simulación reiniciado (ya corrió hasta t_end_years).

        Raises
        ------
        FileNotFoundError
            Si el dump file no existe en datadir.
        ValueError
            Si t_end_years es menor o igual al tiempo actual del dump.

        Warnings
        --------
        ¡Los dump files son objetos pickle! Solo leer dumps propios o de
        fuentes confiables (tripodpy docs).
        """
        import math
        from tripodpy import readdump

        dump_path = os.path.join(datadir, dump_filename)
        if not os.path.isfile(dump_path):
            raise FileNotFoundError(
                f"Dump file no encontrado: {dump_path}\n"
                f"Asegúrate de que la simulación haya corrido al menos un snapshot."
            )

        print(f"\n{'─'*60}")
        print(f"[RESTART] Cargando dump: {dump_path}")
        sim = readdump(dump_path)

        t_now_yr  = float(sim.t) / c.year
        print(f"  → t actual:   {t_now_yr:.3e} yr")
        print(f"  → t objetivo: {t_end_years:.3e} yr")

        if t_now_yr >= t_end_years * 0.9999:
            print(f"  → Simulación ya alcanzó t_objetivo ({t_now_yr:.3e} yr >= {t_end_years:.3e} yr). Nada que correr.")
            print(f"{'─'*60}")
            return sim

        # ── Reconstruir la lista de snapshots ────────────────────────────
        # El dump guarda el array COMPLETO de snapshots originales (incluyendo
        # futuros no completados). Hay que filtrar y extender correctamente.
        t_now_cgs = float(sim.t)

        # Snapshots del plan original que AÚN NO se han completado
        existing_future = sim.t.snapshots[sim.t.snapshots > t_now_cgs]

        t_end_cgs = t_end_years * c.year

        if len(existing_future) > 0 and existing_future[-1] >= t_end_cgs * 0.999:
            # El plan original ya cubre t_end_years — solo usar los pendientes
            sim.t.snapshots = existing_future
            print(f"  → Usando {len(existing_future)} snapshots pendientes "
                  f"del plan original (hasta {existing_future[-1]/c.year:.2e} yr)")
        else:
            # Hay que extender más allá del plan original
            t_max_cgs = existing_future[-1] if len(existing_future) > 0 else t_now_cgs
            t_max_yr  = t_max_cgs / c.year
            t_extra = np.logspace(
                math.log10(t_max_yr),
                math.log10(t_end_years),
                num=num_extra_snapshots + 1,
            )[1:] * c.year   # excluir t_max (ya en el array)
            
            # Prevenir puntos numéricamente idénticos
            t_extra = t_extra[t_extra > t_max_cgs * 1.000001]
            
            sim.t.snapshots = np.concatenate([existing_future, t_extra])
            print(f"  → {len(existing_future)} snapshots pendientes "
                  f"+ {len(t_extra)} nuevos hasta {t_end_years:.1e} yr")

        n_snaps = len(sim.t.snapshots)
        print(f"  → Total snapshots a correr: {n_snaps}")
        print(f"{'─'*60}")

        if n_snaps == 0:
            print("  [!] Sin snapshots pendientes — simulación ya completada.")
            return sim

        # IMPORTANTE: Desactivar overwrite para que simframe continúe 
        # desde el último dataXXXX.hdf5 existente en lugar de empezar en 0000.
        sim.writer.overwrite = False
        sim.run()
        return sim



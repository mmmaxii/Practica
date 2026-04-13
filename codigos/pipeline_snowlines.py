import numpy as np
import dustpy.constants as c
from tripodpy import Simulation

class WaterworldPipeline:
    """
    Pipeline maestro para simular evolución del disco, snowlines (H2O, CO, CO2)
    y acumulación de polvo (Pebble Snow) usando tripodpy, con integración
    a largo plazo sin planetas forzados.
    """
    def __init__(self, datadir="output_test_pipeline_v2"):
        self.sim = Simulation()
        self.datadir = datadir
        
        # Parámetros por default (pueden ser modificados externamente)
        self.M_star_ini = 1.0 * c.M_sun
        self.R_star_ini = 2.0 * c.R_sun  # Radio estelar joven inflado
        self.T_star_ini = 5778.0  # Kelvin (main sequence)

        # ── Especies activas ─────────────────────────────────────────────────
        # Lista de especies de chem.txt que se inyectarán como componentes.
        # Modificar ANTES de llamar a add_volatile_components().
        #
        # Ejemplos:
        #   pipeline.active_species = ["H2O", "CO2", "CO"]          # default
        #   pipeline.active_species = ["H2O", "CO2", "CO", "N2"]    # + nitrógeno
        #   pipeline.active_species = ["H2O"]                        # solo agua
        #
        # Todas las especies deben estar definidas en chem.txt.
        # Las especies con nu_des <= 0 se añaden como trazadores solo-gas.
        # Las especies con nu_des >  0 se añaden como componentes híbridos
        # (condensación/evaporación activa).
        self.active_species = ["H2O", "CO2", "CO"]

        # ── Velocidades de fragmentación para el hielo ──────────────────────
        # Para cada especie condensable, define el umbral de temperatura [K]
        # y la v_frag [cm/s] aplicada cuando T < T_sub.
        # setup_physics() usa esto junto con active_species para construir
        # el updater dinámicamente: solo las especies en active_species
        # aportan un escalón al perfil v_frag(T).
        #
        # Silicatos (baseline, siempre activo): v_frag = 100 cm/s = 1 m/s
        #   Referencia: Birnstiel et al. (2012), A&A 539, A148
        # H2O ice: v_frag = 1000 cm/s = 10 m/s
        #   Referencia: Gundlach & Blum (2015), ApJ 798, 34
        # CO2 ice: v_frag = 500 cm/s = 5 m/s
        #   Referencia: Gundlach et al. (2018), MNRAS 479, 1273
        # CO ice:  v_frag = 300 cm/s = 3 m/s
        #   Referencia: Dominik & Tielens (1997), ApJ 480, 647
        self.vfrag_params = {
            #  especie : (T_sub [K], v_frag [cm/s])
            "H2O" : (150.0, 1000.0),
            "CO2" : ( 70.0,  500.0),
            "CO"  : ( 25.0,  300.0),
        }
        self.vfrag_silicates = 100.0   # baseline siempre activo [cm/s]
        
    def setup_grid(self, rmin=1.0*c.au, rmax=1000.0*c.au, Nr=300):
        """
        Configura la grilla espacial.
        Posteriormente podemos agregar 'grid refinement' de dustpylib aquí si las
        snowlines crean gradientes demasiado fuertes numéricamente.
        """
        print(f"Configurando Grilla: {Nr} celdas entre {rmin/c.au} y {rmax/c.au} AU")
        self.sim.ini.grid.rmin = rmin
        self.sim.ini.grid.rmax = rmax
        self.sim.ini.grid.Nr = Nr
        self.sim.ini.gas.SigmaRc = 30.0 * c.au  # Radio característico (Lynden-Bell & Pringle 1974, R_c ~30 AU típico)
        
    def setup_star(self):
        """
        Inicializa parámetros termodinámicos de la estrella central.
        """
        self.sim.ini.star.M = self.M_star_ini
        self.sim.ini.star.R = self.R_star_ini
        self.sim.ini.star.T = self.T_star_ini
        
    def initialize_simulation(self):
        """
        Inicializa trípodpy. DEBE ser llamado antes de agregar componentes ad-hoc.
        """
        print("Inicializando la simulación base...")
        # Configurar limites físicos y velocidades referenciales
        # sim.ini.dust.vfrag (podemos definir esto de manera global o component by component)
        
        self.sim.initialize()
        
    def add_volatile_components(self):
        """
        Inyecta los trazadores térmicos de volátiles a partir de chem.txt

        Sigue el patrón de referencia de tripodpy:
          - Componentes solo-gas (nu_des <= 0): gas_active=True, dust_active=False
          - Componentes híbridos (nu_des > 0):  gas_active=True, dust_active=True
              → inicializar gas.pars.nu, gas.pars.Tsub, gas.pars.mu y dust.pars.rhos
                explícitamente después de addcomponent (patrón de referencia)
          - Componente de fondo silicatos: dust_active=True, gas_active=False
          - Llamar dust.rhos.update() después de agregar todos los componentes
            para que tripodpy recalcule la densidad interna del polvo correctamente
        """
        print("Agregando sub-componentes (volátiles y refractarios) desde chem.txt...")
        
        data = np.genfromtxt(
            'chem.txt',
            dtype=None,
            names=True,
            encoding='utf-8',
            usecols=(0, 1, 2, 3, 4, 5),
            comments='#'
        )

        # ── Gas residual: se irá descontando a medida que se asignen especies ──
        Sig_residual = self.sim.gas.Sigma.copy()
        frac_h = 0.9118
        self.sim.gas.mu = np.ones_like(self.sim.gas.mu) * (
            frac_h * c.m_p + (1 - frac_h) * 4 * c.m_p
        ) / (0.5 * frac_h + (1 - frac_h))

        mHe = (2 * self.sim.gas.mu/c.m_p - 4) / (4 - self.sim.gas.mu/c.m_p)

        print(f"  → Especies activas: {self.active_species}")

        for element in data:
            spec_name = element['Species']
            if spec_name not in self.active_species:
                continue

            # Fracción de masa relativa a la masa atómica media del disco
            mass_frac = element['Abundance'] * element['mu'] / (frac_h + 4*(1-frac_h))

            if element["nu_des"] <= 0:
                # ── Componente solo-gas (trazador) ───────────────────────────
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma * mass_frac,
                    element['mu'] * c.m_p,
                    dust_active=False,
                    gas_active=True,
                    
                )
                # NO se descuenta aquí: la sustracción única está al final del loop

            else:
                # ── Componente híbrido (condensación/evaporación) ─────────────
                # Truncar la distribución inicial donde el polvo ya está en el floor
                factor = np.where(
                    self.sim.dust.Sigma[:, 0] > self.sim.dust.SigmaFloor[:, 0],
                    1.0, 0.0
                )
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma * mass_frac * factor,
                    element['mu'] * c.m_p,
                    dust_active=True,
                    gas_active=True,
                    rhos=1.,
                    dust_value=self.sim.dust.SigmaFloor.copy(),
                )

                # AQUI debe ser bueno agregar los tracer para trackear la evolucion de los sigma
                # de los componentes.


                # Parámetros de condensación/evaporación — deben setearse
                # DESPUÉS de addcomponent (patrón de referencia de tripodpy)
                comp = getattr(self.sim.components, spec_name)
                comp.gas.pars.nu    = element['nu_des']      # frecuencia de desorción [Hz]
                comp.gas.pars.Tsub  = element['T_bind']      # energía de enlace [K]
                comp.gas.pars.mu    = element['mu'] * c.m_p  # masa molecular [g] — CGS
                comp.dust.pars.rhos = 1.0                    # bulk density polvo [g/cm³]

            # ── Sustracción única (evita doble descuento del branch gas-only) ──
            Sig_residual -= self.sim.gas.Sigma * mass_frac

        # ── Componente de fondo: silicatos refractarios ───────────────────────
        self.sim.addcomponent(
            "silicates",
            self.sim.gas.SigmaFloor,
            1.,
            dust_active=True,
            gas_active=False,
            rhos=3.5,
            dust_value=self.sim.dust.Sigma.copy(),
        )

        # ── Default: gas residual (H₂) ────────────────────────────────────────
        self.sim.components.Default.gas.Sigma   = Sig_residual * 2  # TODO: revisar factor 2
        self.sim.components.Default.gas.pars.mu = 2 * c.m_p

        # ── Actualizar densidad interna del polvo DESPUÉS de añadir todos ─────
        # Necesario porque addcomponent con rhos modificado no actualiza el campo
        # dust.rhos automáticamente (patrón de referencia: sim.dust.rhos.update())
        self.sim.dust.rhos.update()
        self.sim.update()
        
    def setup_physics(self):
        """
        Construye el updater de v_frag(T) dinámicamente a partir de
        self.active_species y self.vfrag_params.

        Solo las especies en active_species Y en vfrag_params contribuyen
        un escalón al perfil. Las demás regiones de temperatura usan el
        baseline de silicatos (self.vfrag_silicates = 100 cm/s).

        Ejemplos:
          active_species=["H2O"]           → escalón solo en T < 150 K (1000 cm/s)
          active_species=["H2O","CO2"]     → escalón en T<150 (1000) y T<70 (500)
          active_species=["H2O","CO2","CO"]→ perfil completo de 3 escalónes
        """
        print("Configurando física dinámica del hielo (v_frag updater)...")

        # Filtrar: solo species que tengan entrada en vfrag_params
        ice_species = [
            sp for sp in self.active_species
            if sp in self.vfrag_params
        ]

        # Ordenar por T_sub DESCENDENTE: aplica primero la especie de mayor T
        # (H2O → CO2 → CO), de modo que cada np.where subsiguiente solo
        # sobreescribe las regiones de menor temperatura.
        ice_species_sorted = sorted(
            ice_species,
            key=lambda sp: self.vfrag_params[sp][0],
            reverse=True
        )

        # Registrar qué escalónes se activaron
        for sp in ice_species_sorted:
            Tsub, vf_ice = self.vfrag_params[sp]
            print(f"  → {sp}: T < {Tsub:.0f} K → v_frag = {vf_ice:.0f} cm/s "
                  f"({vf_ice/100:.0f} m/s)")
        if not ice_species_sorted:
            print("  → Sin especies heladas activas: v_frag = "
                  f"{self.vfrag_silicates:.0f} cm/s (silicatos) en todo el disco")

        # Capturar los parámetros en el closure
        _ice  = [(self.vfrag_params[sp][0], self.vfrag_params[sp][1])
                 for sp in ice_species_sorted]
        _base = self.vfrag_silicates

        def v_frag_variable(sim):
            """
            Perfil v_frag(T) construido dinámicamente.
            Baseline: silicatos = _base cm/s.
            Cada especie helada activa añade un escalón: T < T_sub → v_ice.

            Algoritmo (de mayor a menor T_sub):
              vf = _base
              vf = where(T < T_H2O, v_H2O, vf)   → sobreescribe zona fría
              vf = where(T < T_CO2, v_CO2, vf)    → sobreescribe zona más fría
              vf = where(T < T_CO,  v_CO,  vf)    → sobreescribe la zona más fría
            """
            T  = sim.gas.T
            vf = np.full_like(T, _base)            # baseline silicatos
            for T_sub, v_ice in _ice:              # de mayor a menor T_sub
                vf = np.where(T < T_sub, v_ice, vf)
            return vf

        self.sim.dust.v.frag.updater.updater = v_frag_variable
        self.sim.dust.v.frag.update()

    def add_snowline_fields(self):
        """
        Agrega Fields de posición de snowline al grid para cada volátil activo.

        Sigue el patrón de la documentación de tripodpy:

            sim.grid.addfield("rsnow", 0., description="Snowline Location [cm]")

            def rsnow(sim):
                isnow = np.argmax(sim.gas.T < 150.)
                return sim.grid.ri[isnow]

            sim.grid.rsnow.updater.updater = rsnow

        Se crean tres fields independientes: rsnow_H2O, rsnow_CO2, rsnow_CO.
        Cada uno lee su Tsub directamente del componente correspondiente
        (sim.components.X.gas.pars.Tsub), garantizando consistencia fisica.

        Estos Fields tienen save=True por defecto (es el default de addfield),
        asi que se guardan en cada HDF5 automaticamente con shape ()
        -> despues de read.all() tienen shape (Nt,).

        DEBE llamarse despues de add_volatile_components() y setup_physics().
        """
        print("Agregando fields de posicion de snowline al grid...")

        # Especies con componente activo en la simulacion
        snowline_species = {
            "H2O": 150.0,
            "CO2": 70.0,
            "CO":  25.0,
        }

        def make_rsnow_updater(spec_name, T_sub_K):
            """
            Closure que captura el nombre de la especie y su temperatura de sublimacion.

            IMPORTANTE: comp.gas.pars.Tsub almacena T_bind (energia de enlace en K:
            5800 para H2O, 2700 para CO2, 1180 para CO), NO la temperatura de
            sublimacion fisica. Por eso usamos directamente T_sub_K (150/70/25 K),
            que es lo que define la posicion del snowline.

            La temperatura de sublimacion real se calcula internamente por tripodpy
            usando la ecuacion de Clausius-Clapeyron:
                Tsub ~ T_bind / log(nu * sqrt(2*pi*m*kB*T_bind) / (n_gas * kB))
            pero para propositos de diagnostico, los valores estandar son suficientes.
            """
            def _rsnow_updater(sim):
                # Buscar el primer radio donde T < T_sub (de adentro hacia afuera)
                # Si T[0] < T_sub (disco muy frio): argmax retorna 0 -> ri[0]
                # Si T nunca < T_sub (disco muy caliente): snowline fuera del grid
                T = sim.gas.T
                cold_mask = T < T_sub_K

                if not cold_mask.any():
                    # Disco mas caliente que T_sub en todo el grid:
                    # el snowline esta MAS ALLA del borde exterior
                    return float(sim.grid.ri[-1])

                # Primer indice donde T < T_sub (desde el interior)
                i_snow = int(np.argmax(cold_mask))
                return float(sim.grid.ri[i_snow])

            return _rsnow_updater

        for spec_name, T_sub_default in snowline_species.items():
            field_name = f"rsnow_{spec_name}"
            self.sim.grid.addfield(
                field_name,
                0.0,
                updater=make_rsnow_updater(spec_name, T_sub_default),
                description=f"Snowline position of {spec_name} [cm]",
                save=True
            )
            # Calcular posicion inicial

            getattr(self.sim.grid, field_name).update()
            r_now = float(getattr(self.sim.grid, field_name)) / c.au
            print(f"  -> rsnow_{spec_name} (Tsub={T_sub_default}K): {r_now:.2f} AU")

    def setup_star_evolution(self):

        print("Configurando evolucion estelar (contraccion del radio)...")

        def R_star_evolving(sim):
            # Valores basados en modelos de pre-secuencia principal para 1 M_sun.
            #
            # R_ini: consistente con self.R_star_ini = 2.0 R_sun definido en __init__.
            #   El valor anterior era 3.0 R_sun → causaba una discontinuidad al
            #   activar este updater (la estrella saltaba de 2 a 3 R_sun en t=0+dt).
            #   Referencia: Baraffe et al. (2015), A&A 577, A42 — track de 1 M_sun
            #   a ~1 Myr de edad: R ≈ 1.8–2.2 R_sun.
            #
            # R_fin: la estrella alcanza la ZAMS (~1 R_sun) en ~30 Myr. Dentro
            #   del rango de simulación (10^5 yr), el cambio es mínimo (~1%).
            #   Se fija en 1.5 R_sun como valor intermedio físicamente razonable
            #   durante la fase T Tauri clásica (Baraffe et al. 2015).
            #
            # t_contract: la escala de tiempo de Kelvin-Helmholtz para 1 M_sun
            #   es t_KH ~ G M² / (R L) ~ 10–30 Myr. Se usa 10 Myr como cota
            #   inferior conservadora (Hayashi 1961; Siess et al. 2000).
            #   El valor anterior era 10,000 años → físicamente incorrecto (5 órdenes
            #   de magnitud más rápido que la contracción real).
            R_ini       = 2.0 * c.R_sun    # t~1 Myr (Baraffe+ 2015, 1 M_sun)
            R_fin       = 1.5 * c.R_sun    # fase T Tauri clásica
            t_contract  = 1.0e7 * c.year   # t_KH ~ 10–30 Myr para 1 M_sun (Hayashi 1961; Siess et al. 2000)

            slope = (R_fin - R_ini) / t_contract
            R = R_ini + slope * sim.t
            return float(np.maximum(R, R_fin))

        self.sim.star.R.updater.updater = R_star_evolving
        self.sim.star.R.update()


    def add_ice_sigma_fields(self):
        """
        Agrega fields de densidad superficial de polvo y gas por componente.

        PROBLEMA con updaters en addfield (versión anterior):
          sim.grid fields con updaters se resuelven fuera del orden correcto
          de la cadena de update de tripodpy. El writer llama sim.update(), pero
          'grid' se actualiza antes de que sim.components tenga los valores del
          timestep actual → los fields siempre capturan el estado INICIAL.

        SOLUCIÓN: diastole en sim.dust.updater.diastole
          El diastole se dispara DESPUÉS de que sim.dust (y sim.components)
          termina de calcular el estado del timestep. Asignamos los valores
          directamente in-place con field[:] = value, que sí actualiza el array
          subyacente del Field de simframe.

        Fields creados (para cada especie en active_species):
          grid/SigmaDust_{sp}  — densidad superficial de polvo [g/cm²]  (Nr,)
          grid/SigmaGas_{sp}   — densidad superficial de gas  [g/cm²]  (Nr,)

        Después de read.all() → shape (Nt, Nr) automáticamente.
        """
        print("Agregando fields de densidad superficial por componente (SigmaDust/SigmaGas)...")

        Nr = len(self.sim.grid.r)

        # ── Registrar fields sin updater (se actualizarán vía diastole) ──────
        for sp_name in self.active_species:
            self.sim.grid.addfield(
                f"SigmaDust_{sp_name}",
                np.zeros(Nr),
                updater=None,
                description=f"Σ_dust del componente {sp_name} [g/cm²]",
                save=True,
            )
            self.sim.grid.addfield(
                f"SigmaGas_{sp_name}",
                np.zeros(Nr),
                updater=None,
                description=f"Σ_gas del componente {sp_name} [g/cm²]",
                save=True,
            )

        # ── Diastole: asigned directamente después de cada update de dust ─────
        # Capturar la lista de especies en el closure (inmutable en este punto)
        _active = list(self.active_species)

        def _update_sigma_comp(sim):
            """
            Escribe los valores actuales de comp.dust.Sigma y comp.gas.Sigma
            en los fields del grid usando asignación in-place ([:]).
            Se llama automáticamente después de cada sim.dust.update().

            Se aplica SigmaFloor (igual que tripodpy internamente) para evitar
            que valores sub-numéricos del solver contaminen los cálculos de
            composición en PebbleAccretion2.
            """
            dust_floor = sim.dust.SigmaFloor.sum(-1)   # (Nr,) — piso del polvo total
            gas_floor  = sim.gas.SigmaFloor             # (Nr,) — piso del gas total

            for sp_name in _active:
                comp = getattr(sim.components, sp_name)

                # Polvo: (Nr, Nbins) → suma sobre bins → clamp al piso
                dust_val = comp.dust.Sigma.sum(-1)
                getattr(sim.grid, f"SigmaDust_{sp_name}")[:] = \
                    np.maximum(dust_val, dust_floor)

                # Gas: (Nr,) → clamp al piso
                gas_val = comp.gas.Sigma
                getattr(sim.grid, f"SigmaGas_{sp_name}")[:] = \
                    np.maximum(gas_val, gas_floor)

        self.sim.dust.updater.diastole = _update_sigma_comp

        # Inicializar con valores del estado actual (antes de run)
        _update_sigma_comp(self.sim)
        for sp_name in self.active_species:
            sig_d = float(getattr(self.sim.grid, f"SigmaDust_{sp_name}").max())
            sig_g = float(getattr(self.sim.grid, f"SigmaGas_{sp_name}").max())
            print(f"  → SigmaDust_{sp_name}: max = {sig_d:.3e} g/cm²")
            print(f"  → SigmaGas_{sp_name}:  max = {sig_g:.3e} g/cm²")


    def run_integration(self, t_end_years=5, num_snapshots=30):
        print(f"Configurando integración temporal hasta: {t_end_years:.1e} años con {num_snapshots} snapshots.")
        self.sim.t.snapshots = np.logspace(3, 5, num=num_snapshots, base=10) * c.year
        self.sim.writer.datadir = self.datadir
        self.sim.writer.overwrite = True
        
        print("Empezando ejecución...")
        # Descomentar para correr en la vida real
        self.sim.run()

if __name__ == "__main__":
    pipeline = WaterworldPipeline("data_post_pipeline/pipeline_v3_Sigma_update")
    pipeline.active_species = ["H2O", "CO2", "CO"]
    pipeline.setup_grid(rmin=1*c.au, rmax=300*c.au, Nr=200)
    pipeline.setup_star()
    pipeline.initialize_simulation()
    pipeline.add_volatile_components()
    pipeline.setup_physics()
    pipeline.setup_star_evolution()    # radio estelar contrae -> snowlines migran
    pipeline.add_snowline_fields()     # rsnow_H2O/CO2/CO en HDF5
    pipeline.add_ice_sigma_fields()    # SigmaIce_H2O/CO2/CO/silicates en HDF5
    pipeline.sim.update()              # Por seguridad


    pipeline.run_integration(t_end_years=5, num_snapshots=50)

    print("\nPipeline completado.")



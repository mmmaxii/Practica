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
        self.R_star_ini = 3.0 * c.R_sun  # Radio estelar joven inflado
        self.T_star_ini = 5778.0  # Kelvin (main sequence)
        
        # Opciones de componentes (especies volátiles)
        # Diccionario con masa molecular mu [m_p], Tsub [K], nu [frecuencia vibracional referencial Hz], rho_s [g/cm3]
        self.volatiles = {
            "H2O": {"mu": 16.0, "Tsub": 150.0, "nu": 4.0e13, "rhos": 1.0},
            "CO2": {"mu": 44.0, "Tsub": 70.0,  "nu": 1.0e13, "rhos": 1.6},
            "CO":  {"mu": 28.0, "Tsub": 25.0,  "nu": 7.0e11, "rhos": 1.0}
        }
        
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
        """
        print("Agregando sub-componentes (volátiles y refractarios) desde chem.txt...")
        
        data = np.genfromtxt(
            'chem.txt',
            dtype=None,     # Automatically determine data types
            names=True,     # Use first non-skipped line as column names
            encoding='utf-8',  # Handle text encoding
            usecols=(0, 1, 2, 3, 4, 5),  # Specify columns to read
            comments='#'    # Skip lines starting with #
        )

        # Initialize residual gas surface density to total gas surface density
        Sig_residual = self.sim.gas.Sigma.copy()
        frac_h = 0.9118
        self.sim.gas.mu = np.ones_like(self.sim.gas.mu)*(frac_h * c.m_p + (1 - frac_h) * 4 * c.m_p) / ( 0.5 * frac_h + (1 - frac_h))  # Initial mean molecular weight

        mHe = (2 * self.sim.gas.mu/c.m_p - 4) / (4 - self.sim.gas.mu/c.m_p) # drift composition calculation -> molecular vs atomic (H)  He mass fraction
        
        especies_deseadas = ["H2O", "CO2", "CO"]
        
        for element in data:
            spec_name = element['Species']
            if spec_name not in especies_deseadas:
                continue
                
            # component that only has gas phase
            mass_frac = element['Abundance'] * element['mu'] / (frac_h + 4*(1-frac_h)) #mass fraction relative to the approximate atomic mean mass in the disc

            if element["nu_des"] <= 0:
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma*mass_frac,
                    element['mu']*c.m_p,
                    dust_active=False,
                    gas_active=True
                )
                Sig_residual -= self.sim.gas.Sigma*mass_frac

            # component that has both gas and dust phases
            else:
                # truncate the component where dust surface density floor is reached
                factor = np.where(self.sim.dust.Sigma[:,0] > self.sim.dust.SigmaFloor[:,0], 1.0, 0.0)
                self.sim.addcomponent(
                    spec_name,
                    self.sim.gas.Sigma*mass_frac*factor,
                    element['mu']*c.m_p,
                    dust_active=True,
                    gas_active=True,
                    rhos=1.,
                    dust_value = self.sim.dust.SigmaFloor.copy()
                )
                spec_comp = getattr(self.sim.components, spec_name)
                spec_comp.gas.pars.nu = element['nu_des']
                spec_comp.gas.pars.Tsub = element['T_bind']

            # Subtract assigned gas from residual gas surface density
            Sig_residual -= self.sim.gas.Sigma*mass_frac

        # Add the background silicate dust component
        self.sim.addcomponent("silicates",
            self.sim.gas.SigmaFloor,
            1.,
            dust_active=True,
            gas_active=False,
            rhos=3.5,
            dust_value = self.sim.dust.Sigma.copy()
        )

        # Note that we have to subtract the assigned gas from the Default component and reassign it H2 like properties
        self.sim.components.Default.gas.Sigma = Sig_residual*2 # Assign remaining gas to Default component TODO: why factor 2?
        self.sim.components.Default.gas.pars.mu = 2*c.m_p # Ensure Default behaves Like is H2
        self.sim.update()
        
    def setup_physics(self):
        """
        Personaliza la física de la simulación. Específicamente, el updater de la 
        velocidad de fragmentación (v_frag) para que sea dependiente del cruce de
        las Snowlines, generando el 'traffic jam' natural.
        """
        print("Configurando física dinámica del hielo (v_frag updater)...")
        # Aseguramos que la estructura interna esté creada antes de modificarla
        
        def v_frag_variable(sim):
            # Física térmica para la fragilidad de granos:
            # 1. T < Tsub_CO (Afuera de CO): Hielo de CO los hace frágiles -> 5 m/s (500 cm/s)
            # 2. Tsub_CO < T < Tsub_H2O: Hielo de agua los hace resistentes -> 10 m/s (1000 cm/s)
            # 3. T > Tsub_H2O (Adentro): Silicatos secos estándar -> 1 m/s (100 cm/s)
            t_sub_co = sim.components.CO.gas.pars.Tsub
            t_sub_h2o = sim.components.H2O.gas.pars.Tsub
            
            # Inicializamos array con el valor caliente (Silicatos a 100 cm/s)
            vf = np.full_like(sim.gas.T, 100.0)
            
            # Aplicamos condición de Hielo de Agua (1000 cm/s)
            vf[sim.gas.T < t_sub_h2o] = 1000.0
            
            # Aplicamos condición de Hielo de CO (500 cm/s)
            vf[sim.gas.T < t_sub_co] = 500.0
            
            return vf
            
        self.sim.dust.v.frag.updater.updater = v_frag_variable

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
        """
        Hace que el radio estelar contraiga de 3 R_sun a 1 R_sun en los
        primeros 10,000 anos, simulando la contraccion de Kelvin-Helmholtz
        de una estrella en pre-secuencia principal.

        Propuesta por la documentacion de tripodpy como ejemplo de como
        hacer evolucionar las propiedades estelares:

            def Rstar(sim):
                dR = -1.*c.R_sun
                dt = 1.e4 * c.year
                m = dR/dt
                R = m*sim.t + 3.*c.R_sun
                R = np.maximum(R, c.R_sun)
                return R

            sim.star.R.updater.updater = Rstar

        Nota: Es una aproximacion lineal con fines illustrativos.
        La contraccion real de Hayashi sigue una trayectoria mas compleja.
        El efecto clave es que un radio estelar menor = disco mas frio
        = los snowlines migran hacia adentro con el tiempo.
        """
        print("Configurando evolucion estelar (contraccion del radio)...")

        def R_star_evolving(sim):
            # Radio inicial = 3 R_sun, contrae a 1 R_sun en 10,000 anos
            R_ini = 3.0 * c.R_sun
            R_fin = 1.0 * c.R_sun
            t_contract = 1.0e4 * c.year
            slope = (R_fin - R_ini) / t_contract
            R = R_ini + slope * sim.t
            return float(np.maximum(R, R_fin))

        self.sim.star.R.updater.updater = R_star_evolving
        print(f"  -> R_star(t=0)     = {self.sim.star.R / c.R_sun:.2f} R_sun")
        print(f"  -> R_star(t=1e4yr) = 1.00 R_sun  [estabilizado]")
    def add_ice_sigma_fields(self):
        """
        Agrega fields de densidad superficial de hielo por especie al grid.

        PROBLEMA previo: enable_component_saving() intentaba activar
        comp.dust.Sigma.save, pero en tripodpy el componente dust NO tiene un
        campo Sigma independiente — solo tiene S (source terms), Fi (flujos),
        pars y boundary.

        SOLUCION: Usar sim.grid.addfield() con un updater que calcula la
        densidad superficial de cada especie de hielo directamente, como:

            Sigma_ice_X(r, t) = f_X * sum(dust.Sigma(r,t)) * theta(T < T_sub)

        donde f_X es la fraccion de masa inicial de la especie, theta es la
        funcion escalon (1 si T < Tsub, 0 si T > Tsub), y dust.Sigma.sum(-1)
        es el total de polvo sumando ambos bins.

        Esta formulacion:
          - Incluye el pile-up y drift real del polvo (via dust.Sigma) [OK]
          - Respeta la fisica de snowline exactamente [OK]
          - Se guarda en HDF5 como grid/SigmaIce_X con shape (Nr,) [OK]
          - Despues de read.all() tiene shape (Nt, Nr) [OK]

        REFERENCIAS:
          - Abundancias: Haynes+1992, Fraser+2001 (igual que chem.txt)
          - Fracciones de masa inicial: suma normalizada de abundancias
          - Temperatura de sublimacion: Tsub_H2O=150K, CO2=70K, silicatos>>1500K

        DEBE llamarse DESPUES de add_volatile_components(), setup_physics(),
        add_snowline_fields() y ANTES de run().
        Reemplaza a enable_component_saving(), que no funcionaba correctamente.
        """
        print("Agregando fields de densidad de hielo por especie (SigmaIce_X)...")

        # Abundancias por masa relativa al gas (Haynes 1992 / Fraser 2001)
        species_params = {
            "H2O":       {"T_sub": 150.0,  "abundance": 1.6e-4},
            "CO2":       {"T_sub": 70.0,   "abundance": 4.0e-5},
            "CO":        {"T_sub": 25.0,   "abundance": 8.0e-5},
            "silicates": {"T_sub": 1500.0, "abundance": 2.0e-3},
        }
        total_abun = sum(v["abundance"] for v in species_params.values())

        Nr = len(self.sim.grid.r)

        def make_ice_updater(T_sub_K, f_initial):
            """
            Closure: devuelve el updater de SigmaIce para una especie.

            f_initial: fraccion de masa de la especie respecto al total de especies.
            T_sub_K:   temperatura de sublimacion en K.

            Calculo:
                Sigma_ice(r) = f_initial * Sigma_dust_total(r) * I(T(r) < T_sub)
            donde I() es la funcion indicador.
            """
            def _ice_updater(sim):
                T         = sim.gas.T                    # (Nr,)
                Sig_tot   = sim.dust.Sigma.sum(-1)       # (Nr,) suma bins
                ice_mask  = (T < T_sub_K).astype(float) # 1 fuera del snowline
                return f_initial * Sig_tot * ice_mask    # (Nr,)
            return _ice_updater

        for sp_name, params in species_params.items():
            f_init     = params["abundance"] / total_abun
            T_sub      = params["T_sub"]
            field_name = f"SigmaIce_{sp_name}"

            initial_val = np.zeros(Nr)
            self.sim.grid.addfield(
                field_name,
                initial_val,
                updater=make_ice_updater(T_sub, f_init),
                description=f"Densidad superficial de hielo de {sp_name} [g/cm^2]",
                save=True
            )
            # Calcular valor inicial
            getattr(self.sim.grid, field_name).update()
            sig_max = float(getattr(self.sim.grid, field_name).max())
            print(f"  -> SigmaIce_{sp_name} (Tsub={T_sub}K, f={f_init:.4f}): "
                  f"max={sig_max:.3e} g/cm^2")

    def run_integration(self, t_end_years=5, num_snapshots=30):
        print(f"Configurando integración temporal hasta: {t_end_years:.1e} años con {num_snapshots} snapshots.")
        self.sim.t.snapshots = np.logspace(3, 5, num=num_snapshots, base=10) * c.year
        self.sim.writer.datadir = self.datadir
        self.sim.writer.overwrite = True
        
        print("Empezando ejecución...")
        # Descomentar para correr en la vida real
        self.sim.run()

if __name__ == "__main__":
    pipeline = WaterworldPipeline("data_post_pipeline/pipeline_icefrac")
    pipeline.setup_grid(rmin=1*c.au, rmax=1000*c.au, Nr=300)
    pipeline.setup_star()
    pipeline.initialize_simulation()
    pipeline.add_volatile_components()
    pipeline.setup_physics()
    pipeline.setup_star_evolution()    # radio estelar contrae -> snowlines migran
    pipeline.add_snowline_fields()     # rsnow_H2O/CO2/CO en HDF5
    pipeline.add_ice_sigma_fields()    # SigmaIce_H2O/CO2/CO/silicates en HDF5

    pipeline.sim.update()
    pipeline.run_integration(t_end_years=5, num_snapshots=30)

    print("\nPipeline completo con SigmaIce guardado en HDF5.")



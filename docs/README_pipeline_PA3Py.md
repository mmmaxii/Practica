# Pipeline de Simulacion de Formacion Planetaria
## tripodpy + PA3Py — Documentacion General

**Proyecto:** Pebble Accretion en Discos Estructurados para la Prediccion de Waterworlds
**Workspace:** `c:\astro\Codigos practica + docs + papers`
**Estado:** Pipeline funcional — Etapa de validacion fisica

---

## 1. Contexto Cientifico

El objetivo del proyecto es entender bajo que condiciones un disco protoplanetario
produce planetas ricos en agua ("Waterworlds") en lugar de planetas rocosos secos.

La idea central es que las **snowlines** (lineas de nieve de H2O, CO2 y CO) actuan
como barreras naturales de presion que:
1. Aumentan la velocidad de fragmentacion del hielo -> granos mas grandes (v_frag efecto por especie)
2. Generan acumulaciones de polvo ("traffic jams") al exterior de la snowline
3. Estas acumulaciones son las fuentes de pebbles que acrecen los embriones planetarios

**No se usan planetas artificiales ni gaps forzados.**
Las estructuras emergen organicamente de la termodinamica.

### Referencias clave
- Ormel (2017): *The Emerging Paradigm of Pebble Accretion* — ecuaciones de regimen
- Drazkowska et al. (2023): *Planet Formation Theory in the Era of ALMA and Kepler*
- Musiolik & Wurm (2019): velocidades de fragmentacion del hielo de agua
- Gundlach & Blum (2015): v_frag H2O ice = 10 m/s
- Birnstiel et al. (2012): v_frag silicatos = 1 m/s

---

## 2. Estructura del Proyecto

```
codigos/
  pipeline_snowlines.py       <- Maestro del disco (tripodpy)
  chem.txt                    <- Base de datos de especies quimicas
  plot_diagnostics.py         <- Visualizacion de resultados del disco
  PA3Py/
    PebbleAccretion3.py       <- Modulo de acrecion de pebbles (post-procesador)
    README.md                 <- Referencia matematica del modulo
    test_pa3.py               <- Script de test rapido
    test_rates.py             <- Script de diagnostico de tasas
  data_post_pipeline/         <- Outputs HDF5 de simulaciones
    pipeline_v3_Sigma_update/ <- Corrida de referencia actual
    rmax300_t5.5_sn50_gradvfrag/
    ...
docs/
  README_pipeline_PA3Py.md    <- Este archivo
  README_plot_diagnostics.md
  ...
```

---

## 3. Etapa 1: Evolucion del Disco — `pipeline_snowlines.py`

### 3.1 Que hace

La clase `WaterworldPipeline` encapsula una simulacion completa de disco protoplanetario
usando `tripodpy` (TwoPopPy sub-module, Birnstiel & Pfeil 2025).

Modela:
- Evolucion viscosa del gas (Pringle 1981)
- Crecimiento y fragmentacion de polvo (TwoPopPy, Birnstiel et al. 2012)
- Condensacion/evaporacion de volatiles multi-especie (H2O, CO2, CO)
- Migracion dinamica de las snowlines segun T(r, t)
- Contraccion estelar (radio estelar decae de 2 Rsun a 1.5 Rsun en 10 Myr)

### 3.2 Configuracion tipica

```python
from codigos.pipeline_snowlines import WaterworldPipeline
import dustpy.constants as c

pipeline = WaterworldPipeline("path_outputdata")
pipeline.active_species = ["H2O", "CO2", "CO"]
pipeline.setup_grid(rmin=1*c.au, rmax=300*c.au, Nr=200)
pipeline.setup_star()
pipeline.initialize_simulation()
pipeline.add_volatile_components()   # carga chem.txt
pipeline.setup_physics()             # crea v_frag(T) dinamico
pipeline.setup_star_evolution()      # radio estelar contrae
pipeline.add_snowline_fields()       # guarda rsnow_H2O/CO2/CO en HDF5
pipeline.add_ice_sigma_fields()      # guarda SigmaDust/Gas por especie en HDF5
pipeline.sim.update()
pipeline.run_integration(t_end_years=1e5, num_snapshots=50)
```

### 3.3 Parametros de velocidad de fragmentacion

Definidos en `WaterworldPipeline.vfrag_params`:

| Especie   | T_sub [K] | v_frag [cm/s] | v_frag [m/s] | Referencia              |
|-----------|-----------|---------------|--------------|-------------------------|
| Silicatos | —         | 100           | 1            | Birnstiel et al. (2012) |
| H2O       | 150       | 1000          | 10           | Gundlach & Blum (2015)  |
| CO2       | 70        | 500           | 5            | Gundlach et al. (2018)  |
| CO        | 25        | 300           | 3            | Dominik & Tielens (1997)|

El perfil `v_frag(T)` es una funcion escalonada construida dinamicamente por
`setup_physics()`. Solo las especies en `active_species` contribuyen al perfil.

### 3.4 Componentes quimicas (chem.txt)

Cada especie en `chem.txt` tiene:
- `Abundance`: fraccion numerica respecto al H2 (e.g., H2O = 1.6e-4)
- `nu_des`: frecuencia de desorcion [Hz] (si <= 0, es solo gas, sin polvo)
- `T_bind`: energia de enlace en K (se usa en ecuacion de Clausius-Clapeyron)
- `mu`: masa molecular [uma]

Especies activas por defecto: **H2O, CO2, CO** (+ silicatos como fondo).

### 3.5 Campos guardados en HDF5

Por cada snapshot HDF5 se guardan:

| Campo                  | Shape       | Descripcion                          |
|------------------------|-------------|--------------------------------------|
| `gas/Sigma`            | (Nr,)       | Densidad superficial de gas [g/cm2]  |
| `gas/T`                | (Nr,)       | Temperatura [K]                      |
| `gas/cs`               | (Nr,)       | Velocidad del sonido [cm/s]          |
| `gas/eta`              | (Nr,)       | Parametro de gradiente de presion    |
| `gas/nu`               | (Nr,)       | Viscosidad cinematica [cm2/s]        |
| `dust/Sigma`           | (Nr, 2)     | Sigma polvo por bin (pequeno/grande) |
| `dust/v/rad`           | (Nr, 5)     | Velocidad radial del polvo [cm/s]    |
| `dust/St`              | (Nr, 5)     | Numero de Stokes                     |
| `dust/s/max`           | (Nr,)       | Tamano maximo de grano [cm]          |
| `grid/r`               | (Nr,)       | Grid radial [cm]                     |
| `grid/rsnow_H2O`       | scalar      | Posicion snowline H2O [cm]           |
| `grid/rsnow_CO2`       | scalar      | Posicion snowline CO2 [cm]           |
| `grid/rsnow_CO`        | scalar      | Posicion snowline CO [cm]            |
| `grid/SigmaDust_H2O`   | (Nr,)       | Sigma polvo del comp. H2O [g/cm2]    |
| `grid/SigmaGas_H2O`    | (Nr,)       | Sigma gas del comp. H2O [g/cm2]      |
| `grid/SigmaDust_CO2`   | (Nr,)       | Sigma polvo del comp. CO2 [g/cm2]    |
| ... (idem CO, silicates)| ...        | ...                                  |
| `t`                    | scalar      | Tiempo del snapshot [s]              |

---

## 4. Etapa 2: Acrecion de Pebbles — `PA3Py/PebbleAccretion3.py`

### 4.1 Que hace

`PebbleAccretionModule3` es un **post-procesador** que lee los HDF5 generados
por `tripodpy` y simula el crecimiento de embriones planetarios bajo el flujo
de pebbles del disco.

No modifica el disco — solo consume el flujo radial de pebbles que el disco produce.

### 4.2 Fisica implementada

Todas las ecuaciones siguen Ormel (2017) y Drazkowska et al. (2023):

**Masa de inicio de acrecion de pebbles (onset):**
```
M_PA_onset = St * eta^3 * M_star       (Ormel 2017, Eq 7.11)
```
Si M < M_onset: regimen Safronov balistico (acrecion de planetesimales).

**Regimen Headwind (M < M_hw/sh):**
```
Mdot_2D_hw = sqrt(8 G M t_stop v_hw) * Sigma_peb   (Ormel 2017, Eq 7.13)
```

**Regimen Shear (M > M_hw/sh):**
```
Mdot_2D_sh = 2 R_H^2 * Omega_K * St^(2/3) * Sigma_peb  (Drazkowska 2023, Eq 5)
```

**Masa de transicion Headwind/Shear:**
```
M_hw/sh = v_hw^3 / (8 G Omega_K St)    (Ormel 2017, Eq 7.9)
```

**Transicion 2D -> 3D (turbulencia):**
```
Mdot_eff = Mdot_2D * b_col / (b_col + H_peb * sqrt(8/pi))  (Ormel 2017, Eq 7.24)
```

**Masa de aislamiento (isolation mass):**
```
M_iso = 25 ME * (H/r / 0.05)^3 * (M_star/Msun)   (Drazkowska 2023, Eq 6)
```

**Capa de pebbles (scale height):**
```
H_peb = H_gas * sqrt(alpha / (alpha + St))   (Ormel 2017, Eq 7.25)
```

### 4.3 Uso basico

```python
from PA3Py.PebbleAccretion3 import PebbleAccretionModule3

# Cargar datos del disco
pam = PebbleAccretionModule3.from_datadir(
    "path_outputdata",
    M_star=1.0,     # masa estelar en Msun
    t_min_yr=0.0    # ignorar snapshots antes de este tiempo
)

# Correr crecimiento de embriones
EMBRYOS = [1.0, 2.0, 3.0, 5.0, 10.0, 20.0]   # posiciones en AU
results = pam.run_growth(EMBRYOS, M0_g=pam.M_EARTH * 1e-3)  # semilla ~0.001 ME

# Ver tabla de resultados
pam.summary(results)
```

### 4.4 Formato de resultados

`run_growth()` devuelve un dict `{r_au: array}` donde cada array tiene shape `(Nt, 7)`:

| Col | Campo       | Unidades |
|-----|-------------|----------|
| 0   | `t`         | s        |
| 1   | `M_core`    | g        |
| 2   | `M_H2O`     | g        |
| 3   | `M_CO2`     | g        |
| 4   | `M_sil`     | g        |
| 5   | `r_snow_H2O`| AU       |
| 6   | `M_iso`     | g        |

### 4.5 Como se calcula la composicion quimica

Si el HDF5 contiene `grid/SigmaDust_{sp}` (generado por `add_ice_sigma_fields()`):
- Se interpola la fraccion de cada especie localmente en la posicion del embrion
- Se normaliza sobre el total de polvo en esa celda

Si NO existen esos campos (simulaciones del pipeline antiguas antiguas):
- Se usa un modelo de snowline binario: exterior a r_snow -> condensa, interior -> evapora
- Abundancias tomadas de `_ABUNDANCES` (H2O: 1.6e-4, CO2: 4e-5, CO: 8e-5)

### 4.7 Salida de `summary()`

```
--------------------------------------------------------------------------------
  r [AU]    M [ME]  M_iso [ME]    H2O%    CO2%    Sil%  Tipo
--------------------------------------------------------------------------------
    2.24    12.058       11.98    20.3     0.0    79.7  [W] Waterworld
    3.64     0.049       17.45    16.5     0.0    81.5  [W] Waterworld
    6.75     0.001       27.50     2.7     0.0    13.3  [R] Rocoso
```

`M_iso` es la masa de aislamiento final al ultimo snapshot — indica cuando el
embrion habria alcanzado el limite de acrecion de pebbles segun el estado del disco.

---

## 5. Visualizacion — `plot_diagnostics.py`

La clase `SnowlineDiagnostics` lee los HDF5 y produce:

| Metodo                | Plot                                               |
|-----------------------|----------------------------------------------------|
| `plot_hovmoller()`    | Diagrama espacio-tiempo de epsilon / Sigma_dust    |
| `plot_hovmoller_rt()` | Multi-panel R(t): Sigma_gas, T, a_max, St, eps     |
| `plot_size_distribution()` | Sigma(a, r) en un snapshot                  |
| `plot_pebble_flux()`  | Flujo de masa Mdot_peb(r, t)                       |
| `plot_profiles()`     | 4-panel: eta, St, a_max, Sigma en un snapshot      |
| `plot_hovmoller_comp()` | Hovmoller por especie quimica                    |

Todas las snowlines se sobreimprimen automaticamente en los plots.

---



## 6. Entorno de Trabajo

**Python requerido:** 3.11.x (no 3.14 — da errores al compilar extensiones C de tripodpy)

**Dependencias:**
- `tripodpy` (instalado como paquete en el venv desde el repo de GitHub)
- `dustpy` (dependencia de tripodpy)
- `numpy`, `matplotlib`, `h5py`

**Activar entorno:**
```powershell
# En la raiz del proyecto
.\env_tripod\Scripts\activate
# o segun donde este el venv
```

**Instalar tripodpy (una sola vez):**
```bash
git clone https://github.com/stammler/tripodpy
cd tripodpy
pip install -e .
```

---

## 8. Limitaciones Actuales del Modelo

| Limitacion              | Impacto                                                    |
|-------------------------|------------------------------------------------------------|
| Embriones fijos en r    | Sin migracion Tipo I — composicion subestimada             |
| Sin acrecion de gas     | No hay casos Super-Tierra/Neptuno — solo cores rocosos     |
| Sin foto-evaporacion    | El disco no se dispersa — M_iso no cae                     |
| 1 estrella fija         | Sin variacion de M_star en la misma corrida                |
| t_end corto (~10^5 yr)  | No llega a 1-5 Myr donde la fisica es mas rica             |

---

## 9. Proximos Pasos (Plan de Validacion)

Ver `validation_plan.md` en el directorio de la IA para el plan detallado.

Resumen:
1. **Check 1** — Flujo de pebbles vs. formula analitica de Lambrechts & Johansen (2014)
2. **Check 2** — St, eta, v_hw vs. rangos esperados de la literatura
3. **Check 3** — M_iso(r) vs. valores de Drazkowska et al. (2023)
4. **Check 4** — Escala de tiempo de acrecion: t(1 ME) ~ 0.1-5 Myr
5. **Check 5** — Composicion quimica vs. posicion de snowlines

# Pipeline de Simulaciones PPOLs (`pipeline_snowlines.py`)

Este documento explica el funcionamiento del pipeline principal utilizado para las simulaciones de discos protoplanetarios y líneas de nieve (snowlines), construido sobre `tripodpy` y `dustpylib`.

## Arquitectura del Pipeline

El núcleo del pipeline es la clase `WaterworldPipeline` (ubicada en `pipeline_snowlines.py`). Está diseñada bajo un patrón de "Mixins" para modularizar la física y configuración del disco. Los módulos principales se encuentran en la carpeta `pipeline_methods/`:

1.  **`DiskSetupMixin` (`disk_setup.py`)**: 
    Se encarga de inicializar la simulación base. Define la estrella central (masa, radio, temperatura), la grilla radial logarítmica y los parámetros iniciales del gas (como `alpha` y masa del disco).
    
2.  **`DiskChemistryMixin` (`disk_chemistry.py`)**: 
    Configura la composición del polvo. Inyecta especies químicas activas (por defecto H2O, CO2, CO, más los silicatos base) leyendo sus propiedades desde `chem.txt`. También maneja la distribución inicial de densidad superficial (Σ) para el hielo.

3.  **`SnowlinePhysicsMixin` (`snowline_physics.py`)**: 
    Implementa la física de las velocidades de fragmentación ($v_{\rm frag}$). Define cómo la velocidad de fragmentación cambia en función de la temperatura del disco, detectando dinámicamente las posiciones de las líneas de sublimación (snowlines) para cada especie.

4.  **`PressureBumpsMixin` (`pressure_bumps.py`)**: 
    Maneja la formación de subestructuras (gaps planetarios y trampas de presión) modificando el perfil de viscosidad del gas (`gas.alpha`).

---

## Flujo de Ejecución Típico

Una simulación estándar se inicializa y ejecuta de la siguiente manera:

```python
from pipeline_snowlines import WaterworldPipeline
import dustpy.constants as c

# 1. Instanciar el pipeline con parámetros estelares y de grilla
pipeline = WaterworldPipeline(
    datadir="mi_simulacion/",
    active_species=["H2O", "CO2", "CO"],
    grid_rmin=1 * c.au,
    grid_rmax=300 * c.au,
    Nr=200,
    M_star_Msun=1.0
)

# 2. (Opcional) Agregar un gap planetario
pipeline.setup_gap_duffell(M_planet=30.0 * c.M_earth, a_planet_au=5.0, imprint=True)

# 3. Iniciar la integración temporal
pipeline.run_integration(t_end_years=1e6, num_snapshots=50)
```

Durante el `__init__`, el pipeline automáticamente:
* Configura la grilla y la estrella.
* Inicializa la simulación base en `tripodpy`.
* Agrega los componentes volátiles.
* Configura la física dependiente de la temperatura y las líneas de nieve.
* Realiza un `update()` final para dejar todo listo para integrar.

---

## Modelado de Gaps Planetarios y Subestructuras

La modificación de la estructura del disco para simular planetas se realiza a través de `PressureBumpsMixin`. Este módulo altera el parámetro de viscosidad $\alpha(r)$ asumiendo que el flujo en estado estacionario cumple $\Sigma \propto 1/\alpha$. Al reducir artificialmente $\alpha$ en ciertas zonas, se generan acumulaciones de gas y polvo (trampas de presión).

Los perfiles se calculan de manera estática en $t=0$. Se recomienda usar la opción `imprint=True` para tallar el disco (modificar $\Sigma_{\rm gas}$ y $\Sigma_{\rm dust}$) inmediatamente y evitar flujos viscosos violentos durante los primeros pasos de integración.

### Métodos Disponibles:

1.  **Modelo de Kanagawa et al. (2017)**
    Abre un gap basado en el torque de un planeta. Útil para planetas masivos (tipo Júpiter).
    ```python
    pipeline.setup_gap_kanagawa(M_planet=1.0 * c.M_jup, a_planet_au=5.2, imprint=True)
    ```

2.  **Modelo de Duffell (2020)**
    Es una alternativa empírica que ajusta bien para planetas de menor masa (ej. Súper-Tierras, Neptunos).
    ```python
    pipeline.setup_gap_duffell(M_planet=30.0 * c.M_earth, a_planet_au=5.0, imprint=True)
    ```

3.  **Múltiples Planetas (Duffell Multi)**
    Permite definir un sistema con varios planetas simultáneamente. Multiplica los perfiles de atenuación para cada planeta.
    ```python
    planets = [
        {"M_planet": 30 * c.M_earth, "a_planet_au": 5.0},
        {"M_planet": 30 * c.M_earth, "a_planet_au": 10.0}
    ]
    pipeline.setup_gap_duffell_multi(planets=planets, imprint=True)
    ```

4.  **Perfil Sinusoidal**
    Genera múltiples gaps equidistantes (en escala logarítmica) perturbando $\alpha$ con una función sinusoidal. Útil para estudiar trampas de presión periódicas (anillos) sin invocar planetas específicos.
    ```python
    pipeline.setup_alpha_sinusoidal(amplitude=5.0, n_bumps=5, r_inner_au=5.0, r_outer_au=100.0, imprint=True)
    ```

### Reinicio de Gaps
Si se necesita comparar o limpiar el estado de la viscosidad en un entorno interactivo (Jupyter), se puede utilizar:
```python
pipeline.reset_gap()
```
Esto restaura el valor de $\alpha$ a su perfil constante original.

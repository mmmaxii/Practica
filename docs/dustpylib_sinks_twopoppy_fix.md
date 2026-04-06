# Integración de Sumideros de DustPyLib en TripodPy (TwoPopulations)

Al utilizar los módulos avanzados de formación de planetesimales de `dustpylib` (como `drazkowska2016`, `miller2021` y `schoonenberg2018`) dentro de un entorno basado en `tripodpy` o la arquitectura Two-Populations, es completamente natural que enfrentes el siguiente error de colisión de arreglos (*Broadcasting ValueError*):

`ValueError: operands could not be broadcast together with shapes (150,2) (150,5)`

---

## ¿Por qué ocurre y qué significa físicamente?

La librería `dustpylib` fue pensada para integrarse directo a un `dustpy` estándar, donde la masa del polvo se distribuye en decenas de contenedores o perfiles masivos continuos ($N_m$ mass-bins). `tripodpy`, a favor de computación rápida, abstrae gran parte del polvo a un modelo `TwoPopPy`. 

- La variable maestra de tus masas **`sim.dust.Sigma`**, se representa lógicamente con una estructura de `(150, 2)`:
  1. Columna `[0]`: Superficie de masas de todo el polvo minúsculo y base.
  2. Columna `[1]`: Superficie masiva de todo el polvo grande visiblemente dominante (los *Pebbles* astronómicos).

- Sin embargo, tu parámetro aerodinámico principal, el Stokes Number **`sim.dust.St`** y tu barrera de tamaño diametral **`sim.dust.a`**, mantienen una forma de vector fantasma o tracker constante de `(150, 5)`. Y su llave fundamental acusa lo siguiente:
  - `a[0]`: $a_0$ (Polvo estelar base)
  - `a[1]`: $\text{fudge} \cdot a_1$ 
  - `a[2]`: $a_1$ 
  - `a[3]`: $0.5 \cdot a_{max}$ 
  - `a[4]`: $a_{max}$ (El polvo de frontera, dominado por límites radiculares u opacidades de fragmentación)

Las ecucaciones formulativas de las funciones de gaps obligaban tu matriz de dos entidades (el Sigma) a operar diferencialmente frente la aerodinámica de los 5 objetos de cálculo del campo estocástico.

---

## La Solución: "Indexación Adaptativa (Fancy Numpy Indexing)"

Sería un peligro astronómico gigante aglomerar al primer pilar de la matriz poblacional (`Sigma` de polvo minúsculo) los comportamientos de choque físico para cálculos planetesimales del pilar equivocado (Por ejemplo de un índice `3`). 

La respuesta para no corromper la termodinámica integrativa es blindar la llamada entregando al actualizador **solo las dos columnas maestras del Campo Predictivo de Stokes ($St$) congruentes con la Pista de Masa Observable ($Sigma$)**:
1. Para atar al índice del polvo pequeño de masa `Sigma[:, 0]`, proveemos los componentes aerodinámicos de la semilla $St_{a0}$: **`St[:, 0]`**
2. Para emparejar a los aglomerados masivos y guiar bien nuestro "gap sumidero", amarramos genéticamente `Sigma[:, 1]` con sus contrapartes pebbles nativas $St_{amax}$: **`St[:, 4]`**

El corte lo hacemos instanciando variables internas con doble corchete, adaptando el formato previo al envío numérico.

---

## Snippets Blindados Universales

Cualquiera de tus *Updaters* instalados a `sim.dust.S.ext.updater` deben ser configurados de la siguiente forma infalible para salvar tu notebook:

**1. Parche para Schoonenberg et al. (2018)**
```python
import numpy as np

def schoo_updater(sim):
    # Selección Inteligente de columnas: Las dos representativas del TwoPopPy
    St_adapted = sim.dust.St[:, [0, 4]] 
    
    return schoonenberg2018(sim.grid.OmegaK, sim.dust.rho, sim.gas.rho, sim.dust.Sigma, St_adapted, d2g_crit=1.0, zeta=0.1)
```

**2. Parche para Miller et al. (2021)**
```python
import numpy as np

def miller_updater(sim):
    # Adaptación
    St_adapted = sim.dust.St[:, [0, 4]] 
    
    return miller2021(sim.grid.OmegaK, sim.dust.rho, sim.gas.rho, sim.dust.Sigma, St_adapted, d2g_crit=1.0, n=0.03, zeta=0.1)
```

**3. Parche para Drążkowska et al. (2016)**
```python
import numpy as np

def draz_updater(sim):
    # Adaptación
    St_adapted = sim.dust.St[:, [0, 4]] 
    
    return drazkowska2016(sim.grid.OmegaK, sim.dust.rho, sim.gas.rho, sim.dust.Sigma, St_adapted, p2g_crit=1.0, St_crit=0.01, zeta=0.01)
```

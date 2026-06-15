# Guía de Uso de DustPyLib en DustPy / TripodPy

Este documento resume las tres áreas funcionales principales de la librería `dustpylib` observadas en su código fuente, los parámetros que involucran y las instrucciones precisas para implementar e inyectar esta física en tu objeto `Simulation`.

---

## 1. Formación de Planetesimales (Sinks / Sumideros)

Estas funciones calculan cuánto polvo (pebbles/granos) se transforma en planetesimales, lo cual computacionalmente se traduce en **sacar masa** de la grilla de polvo dinámico en el disco. Estos modelos asignan tasas de decaimiento negativas `S`.

### Modelos Disponibles:
1. **`drazkowska2016`**: 
   * **Dinámica**: Activa la formación de planetesimales si el ratio polvo-gas en el plano medio (`p2g_mid`) supera un límite crítico (`p2g_crit`). Solo toma en cuenta partículas ("pebbles") con un Stokes number mayor al límite `St_crit`.
   * **Parámetros Únicos**: `p2g_crit` (default=1.), `St_crit` (default=0.01), `zeta` (Eficiencia, default=0.01).
2. **`miller2021`**: 
   * **Dinámica**: Ocupa una transición suave matemática (usando tangente hiperbólica `tanh`) basándose en una evaluación sobre el ratio polvo/gas de todo el continuo.
   * **Parámetros Únicos**: `d2g_crit` (Límite crítico, default=1.), `n` (Parámetro de suavizado suave, default=0.03).
3. **`schoonenberg2018`**: 
   * **Dinámica**: Simple y agresiva; si el ratio polvo-gas supera el crítico local, castiga inmediatamente multiplicando por la eficiencia `zeta` y el Stokes Number.

### ¿Cómo implementarlo?
Cualquiera de estos métodos devuelve una matriz 2D (Radios, Masas) que representa la pérdida de masa externa del polvo. Para inyectarlo, debes **actualizar el campo estocástico de fuentes/sinks (`sim.dust.S.ext`)**.

```python
from dustpylib.planetesimals import drazkowska2016

def planetesimal_updater(sim):
    # Recolectamos la física actual del simulador
    return drazkowska2016(
        sim.grid.OmegaK, 
        sim.dust.rho, 
        sim.gas.rho, 
        sim.dust.Sigma, 
        sim.dust.St, 
        p2g_crit=1.0, 
        St_crit=0.01
    )

# Inicializamos y sobreescribimos
sim.initialize()
sim.dust.S.ext.updater.updater = planetesimal_updater
```

---

## 2. Subestructuras / Planetary Gaps

Estas funciones te permiten moldear los surcos y huecos (`Gaps`) que dejaría un planeta masivo mediante una formulación analítica, alterando la fluidez del gas localmente. 

### Modelos Disponibles:
Ambas funciones devuelven el perfil de perturbación de la densidad superficial $\Sigma / \Sigma_{unperturbed}$.
1. **`kanagawa2017`**: Utiliza el tamaño de la caída del gas y sus paredes de forma lineal/analítica escalando muy bien en la física de "Shallow/Deep Gaps".
2. **`duffell2020`**: Utiliza formulaciones basadas en el número de Mach para esculpir el hueco planetario limitando fallas asintóticas en huecos extremadamente pequeños.

* **Parámetros Universales**: 
  - `r`: El grid radial `sim.grid.r`.
  - `a`: Distancia del planeta a la estrella (Semi-major axis).
  - `q`: Ratio de masa Planet/Estrella ($M_p / M_\star$).
  - `h`: Aspect Ratio (Escala de Altura sobre radio $H/r$) evaluada **en la posición del planeta** o usando el perfil continuo global.
  - `alpha0`: El parámetro de viscosidad $\alpha$ inalterado/base del entorno.

### ¿Cómo implementarlo?
La forma correcta de inducir este hueco no es pisando `Sigma` (lo cual sería inestable a difusiones térmicas), sino que **dividiendo tu parámetro de viscosidad $\alpha$ base por este perfil**, así la bajada de turbulencia esculpe y retiene la masa y trampa natural!

```python
from dustpy import constants as c
from dustpylib.substructures.gaps import kanagawa2017

def gap_updater(sim):
    alpha0 = 1e-3  # Alpha original del disco
    a_planet = 30 * c.au # Planeta en 30 AU
    q_planet = 1e-3 # Ej: Masa de 1 Júpiter
    
    # Aspect Ratio = sim.gas.Hp / sim.grid.r
    h = sim.gas.Hp / sim.grid.r  
    
    # Retorna la perturbación (valor ~1 lejos, valor >1 o <1 dentro)
    f = kanagawa2017(sim.grid.r, a_planet, q_planet, h, alpha0)
    
    # Dividimos alpha original por la perturbación de Kanagawa
    return alpha0 / f

sim.initialize()
sim.gas.alpha.updater.updater = gap_updater
```

---

## 3. Backreaction de Polvo a Gas

Clásicamente en DustPy el gas empuja el polvo a migrar (Radial Drift). Sin embargo, cuando se densifica mucho polvo (Bumps o Ice-Lines masivas), **el polvo comienza a repeler y arrastrar al gas consigo**. Esto se denomina *Backreaction*. `dustpylib` cuenta con un wrapper de inyección en 1 sola línea de código para esto.

La librería tiene dos setups (Garate et al. 2019/2020):
1. **Uniforme (`vertical_setup=False`)**: Se asume que el polvo y gas están perfectamente revueltos a cualquier altura `z`.
2. **Estratificado (`vertical_setup=True`)**: Forma más realista. El polvo más pesado se asienta en el plano transversal (Midplane Settling). Modula las velocidades radiales basándose en perfiles gaussianos estratificados utilizando coeficientes A y B por especie. También mitiga la difusividad (`dustDiffusivity_Backreaction`) considerando la obstrucción de colisiones.

### ¿Cómo implementarlo?
No necesitas inyectar *updaters* manuales, sus desarrolladores proveen una cómoda función `setup_backreaction` para re-linkear la física de `sim.dust.v.rad` y asociar los coeficientes ocultos `sim.dust.backreaction.A` y `B`.

```python
from tripodpy import Simulation
from dustpylib.dynamics.backreaction.setup_backreaction import setup_backreaction

sim = Simulation()
sim.initialize()

# Activa el backreaction avanzado automáticamente, pisando los updaters predeterminados de la simulación
setup_backreaction(sim, vertical_setup=True)

sim.run()
```
---
**Nota sobre Créditos:** Utiliza siempre estas fórmulas pre-programadas importándolas como se instruyó, de esta manera aseguraremos máxima rapidez de computación (ya que internamente actúan utilizando vectorización numpy limpia) ahorrando tiempo al no tener que replicar integrales, enmascarados y constantes por nuestra cuenta.

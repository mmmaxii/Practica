# Generación de Pressure Bumps Modificando la Turbulencia (Parámetro Alpha)

## ¿Cuál es el problema original?
Al intentar generar una "trampa de presión" (pressure bump) para acumular polvo, se nos podría ocurrir añadir artificialmente masa al gas, alterando bruscamente `sim.gas.Sigma`. Sin embargo, esto falla en una simulación estándar porque no tiene sustento viscoso continuo. Debido a la evolución física del disco dictada por la ecuación de difusión (gobernada por el parámetro de viscosidad $\alpha$), cualquier abultamiento repentino de densidad es rápidamente difuminado hacia los lados por la fricción viscosa en muy pocos pasos de integración temporal. El sistema vuelve rápidamente a un punto de equilibrio suavizado.

## La Solución Ortodoxa: Modificar el Parámetro Alpha (Dead Zone)
Para generar una trampa de presión realista, duradera y persistente en el tiempo en un disco de acreción viscoso, el método convencional y físicamente fundamentado es **crear una reducción localizada en el parámetro de turbulencia $\alpha$** (lo que simula una zona muerta o "dead zone" de fricción).

Cuando la viscosidad desciende en una franja radicular elegida, al gas "le cuesta fluir" más por ahí. Para seguir respetando la continuidad del flujo másico, dictado por principios de conservación, la masa del gas **se acumulará de manera completamente natural** al encontrarse con ese atasco, y así construirá una pendiente permanente de densidad. 

Esta barrera de gas genera gradientes termodinámicos ($\eta$) sumamente pronunciados que repelerán de forma natural el migrar centrípeto radial del polvo (_radial drift_), haciéndolo quedar capturado.

## Implementación vía Updaters (Heartbeat de Simframe)
No sirve definir `sim.gas.alpha` como un arreglo fijo una sola vez. Cuando se active la revaluación del motor del objeto de simulación (su mecanismo interno de `Updaters` regidos por la arquitectura `Systole/Diastole`), nuestro parche temporal estático de $\alpha$ sería "pisado" y sobreescrito por la función normal que genera variables default para volver a los rangos de inicio. 

Hay que inyectar nuestra matemática como una entidad dinámica permanente en el núcleo. Reemplazamos la orden por defecto del `updater.updater` de `alpha` para que el motor entienda que, en **cada micro-fotograma temporal de la simulación**, la distribución de turbulencia de nuestro disco debe ser dictaminada por este bajón iterativo (_dip_):

```python
import numpy as np
import dustpy.constants as c
from tripodpy import Simulation

# 1. Definimos una curva de campana (Gaussiana) inversa restada de la fricción base
def alpha_dip_updater(sim):
    base_alpha = 1e-3  # Turbulencia estándar
    bump_alpha = 1e-4  # Turbulencia reducida (la zona con bache)
    
    r0_1 = 30 * c.au  ;  width_1 = 5 * c.au
    r0_2 = 100 * c.au ;  width_2 = 10 * c.au
    
    r = sim.grid.r
    
    dip_1 = (base_alpha - bump_alpha) * np.exp(- (r - r0_1)**2 / (2 * width_1**2))
    dip_2 = (base_alpha - bump_alpha) * np.exp(- (r - r0_2)**2 / (2 * width_2**2))
    
    nuevo_alpha = base_alpha - dip_1 - dip_2
    
    return np.maximum(nuevo_alpha, bump_alpha) # Clamping numérico de prevención

# 2. Inicialización obligatoria
sim = Simulation()
sim.ini.grid.Nr = 150
sim.makegrids()
sim.initialize()

# 3. Intervenimos el actualizador directo del parámetro
sim.gas.alpha.updater.updater = alpha_dip_updater

# Forzamos la primera pasada de la nueva física para asentar valores
sim.gas.alpha.update()
sim.update()

# 4. Soltamos el tiempo
sim.writer.datadir = "bumps/pb_alpha"
sim.run()
```

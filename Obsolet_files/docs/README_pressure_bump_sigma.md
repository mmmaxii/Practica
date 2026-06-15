# Generación de Pressure Bumps Forzando la Integración (Variable Sigma)

## El Problema de Sobre-escribir Variables Maestras
Generar una trampa de presión usando el parámetro de turbulencia $\alpha$ es físicamente la forma más idónea. ¿Pero qué sucede si lo que queremos es derechamente esculpir e inyectar permanentemente un muro de densidad inamovible (una alteración en `sim.gas.Sigma` con un perfil Gaussiano directo) obligando al disco a adaptarse a él, a modo de experimento puramente analítico?

Al principio podríamos pensar en aplicar la misma técnica que con `alpha` e imponer una simple instrucción matemática en un *Updater* de la densidad: `sim.gas.Sigma *= multiplier_bump`.

### ¡Peligro de Explosión Numérica!
A diferencia de un parámetro ordinario como la fricción $\alpha$ o un rango limitante (snowline), el tensor **`gas.Sigma` es la mismísima variable maestra de integración condicional del solver matemático.**
Esto quiere decir que `simframe/tripodpy` utiliza métodos iterativos para resolver la ecuación de difusión sumando los flujos a esta matriz. Si usamos un Updater normal para ordenarle cíclicamente la orden de multiplicar la masa  (`Sigma = Sigma * bump`) en la línea de tiempo, esto sucedería en cada sub-loop... una infinidad de veces por año de simulación. Terminaremos operando multiplicaciones factoriales exponenciales y toda simulación implosionaría numéricamente ante la creación constante e ilógica de un universo infinitamente denso (Runtime errors, Not a Number variables, etc).

## La Solución: "Clamping" sobre el Integrador a través de "Hooks" (Ganchos)
Para forzar una escultura inquebrantable y permanentemente estática del gas (creando una muralla inmutable) eludiendo retro-multiplicaciones destructivas, debemos evadir al integrador con sutileza, alterando creativamente los dominios `Fields` y empleando métodos que superan la inferencia del simulador:

1. **Fotografía Intocable**: Guardamos un _snapshot_ perfecto (copia prístina inamovible independiente) de `Sigma` justo al inicializar, exento del dominio temporal. 
2. **Override Total sobre la Matriz Prístina**: Nuestro calculador se basará entonces SIEMPRE en la fotografía intocable. Prevenimos toda re-multiplicación continua destruyendo los deltas que arroje el avance físico del disco, pues le impondremos a la dura nuestra plantilla gaussiana. Pisotearemos (_clamp_) la respuesta de masa integradora al imponer una estampa inamovible.
3. **El Campo Inyector**: Puesto que `Sigma` es variable de evolución, **carece** lógicamente de un mecanismo simple `updater.updater`. Para resolver esto meteremos un campo fantasma (Dummy field) nulo o "Hook" en el abanico principal de evaluación del Gas a través del `addfield`. Asignaremos aquí nuestra instrucción matemática de machacado y lo acomodaremos artificialmente en lo alto (Puesto Cero, `[0]`) de la jerarquía de prioridades evolutiva o `updateorder`.
A partir de tal intromisión, la densidad del disco será obligatoriamente deformada ante los ojos de simulación, la cual calculará sus presiones, $\eta$, y *Stokes* postizos subordinándose maravillosamente a nuestro muro.

```python
import numpy as np
import dustpy.constants as c
from tripodpy import Simulation

sim = Simulation()
sim.ini.grid.Nr = 150
sim.makegrids()
sim.initialize()

# 1. Guardado imperativo originario: la plantilla fotográfica de control maestro
sigma_base = sim.gas.Sigma.copy() 

# 2. Programación del Clamp Opresivo de masa
def forzar_sigma_bump(sim):
    # Definimos amplitudes absolutas agresivas
    A1 = 5.0  ; r0_1 = 30 * c.au  ; width_1 = 5 * c.au
    A2 = 5.0  ; r0_2 = 100 * c.au ; width_2 = 10 * c.au
    
    r = sim.grid.r
    bump1 = A1 * np.exp(- (r - r0_1)**2 / (2 * width_1**2))
    bump2 = A2 * np.exp(- (r - r0_2)**2 / (2 * width_2**2))
    
    # "PISOTEAMOS" el avance lógico del Integrador Viscoso. No mutamos la iteración 
    # actual, mutamos en base de la imagen muerta intocable: ESTO SALVA LA EVOLUCIÓN.
    sim.gas.Sigma = sigma_base * (1 + bump1 + bump2)
    return 0. # Retorno en falso, no guardamos este valor. En la clase group, Dummy.

# 3. Fabricación del Hook Fantasma "imposicion_densidad" insertado al motor de evaluación Simframe:
sim.gas.addfield("imposicion_densidad", 0.)  # Inyección del Campo Falso
sim.gas.imposicion_densidad.updater.updater = forzar_sigma_bump  # Inserción de funciones

# Subida jerárquica obligatoria [Index 0]! Así se garantiza que rija cada etapa micro del loop 
# del gas antes de que los polvos interactúen o consulten gradientes de presión
sim.gas.updateorder.insert(0, "imposicion_densidad") 

# Actualizamos forzadamente nuestro gancho al primer Frame estático
sim.gas.imposicion_densidad.update()
sim.update()

sim.writer.datadir = "bumps/pb_sigma"
sim.run()
```

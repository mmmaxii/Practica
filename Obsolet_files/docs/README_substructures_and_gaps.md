# Guía Exhaustiva: Creación de Gaps y Subestructuras en DustPy / TripodPy

Tras investigar a fondo la arquitectura subyacente del *Standard Model* de `dustpy` y las implementaciones de referencia de los módulos de `dustpylib` (Kanagawa et al. 2017, Duffel 2020, Drążkowska 2016, etc.), es evidente que existen **3 familias metodológicas** para simular substructuras (gaps planetarios, zonas muertas o anillos). 

Cada método tiene su nivel de intrusión, sus ventajas y su modo de lidiar con los solvers implícitos de la arquitectura de la simulación.

---

## Método 1: Modificación del Parámetro Viscoso $\alpha$ (El Método Físico/Elegante)

Este método se aprovecha del flujo hidrodinámico natural. En un disco gaseoso evolutivo dominado por viscosidad, el flujo másico estacionario rige la relación estricta: $\Sigma \times \nu = \text{Constante}$ (donde $\nu = \alpha \cdot c_s \cdot H_p$).

Por tanto, si perturbamos $\alpha$, el motor de `dustpy` alterará la función hidráulica del disco:
- Si **disminuyes** $\alpha$ ("Dead Zone"), el gas se coagulará, generando un bache natural ($\Sigma$ sube), lo cual altera la presión suavemente y atrapa de forma natural el polvo. *¡Este fue el primer intento exitoso!*
- Si **aumentas** $\alpha$ puntualmente, la fricción expulsará el gas local, tallando un natural **Planet Gap** ($\Sigma$ baja).

**Implementación según `dustpylib` (Kanagawa et al. 2017 / Duffell 2020)**:
1. Usan perfiles analíticos ajustados de simulaciones hidrodinámicas 3D (e.g. la función `kanagawa2017` u `duffell2020`) que calculan la distribución de perturbación inversa de densidad respecto de un planeta con masa gravitatoria.
2. Multiplican el inverso de este perfil a `sim.gas.alpha`: `alpha /= kanagawa2017(...)` usando un Actualizador (`updater`).
3. Opcionalmente (para no tener que esperar milenios a que el disco alcance el estado estacionario viscoso) ajustan `Sigma` al estado final desde $t=0$: `s.gas.Sigma[...] /= s.gas.alpha/alpha0`.

> **[!NOTE]**
> Si se prefiere esta vía, notar que actualmente **no está instalado el paquete `dustpylib`** en el entorno local (`ImportError: No module named 'dustpylib'`). Hay que ejecutar `pip install dustpylib` en la terminal activa dentro de `C:\astro\tripodpy\env_tripod\`.

---

## Método 2: Inyección de Sumideros Externos (`S.ext`) (El Método de Acumulación)

Un método nativo para las partículas sólidas y el polvo es alterar su fuente y sumidero (`Sources & Sinks`). En lugar de alterar la hidrodinámica del gas, se retira el polvo a medida que se acumula o alcanza ciertos umbrales, generando un hueco real en el inventario de la distribución de masas.

**Cómo funciona en `dustpylib.planetesimals` (Drążkowska, Schoonenberg, Miller)**:
1. Se evalúa en cada iteración la matriz `sim.dust.S` (Sources). Específicamente el campo de fuerzas externas `sim.dust.S.ext`.
2. Se asigna un `updater` a ese campo que evalúe lógicas físicas; por ejemplo: *"Si localmente la relación Polvo / Gas supera 1, restar masa equivalente a una constante probabilística"* (convertirlo en cuerpos más rígidos "planetesimales", que ya no obedecen el arrastre del gas).
3. Esto inyecta números **negativos** en `S.ext`, lo que remueve masa directamente de `sim.dust.Sigma` de forma consistente durante la matriz de resolución de Euler.
4. Esa masa negativizada se transfiere sumándola a un campo estático explícito (por ejemplo, uno generado localmente mediante `sim.integrator.instructions.append(Instruction(schemes.expl_1_euler, ...))`).

**Uso Sugerido:** Exclusivo para física de acreción (transformación de "Pebbles a Planetesimals") donde ocurre una desaparición genuina de materia del rango dinámicamente acoplado del gas.

---

## Método 3: Imposición Forzosa en Variable de Integración $\Sigma$ (El "Hacking" Computacional)

Esta es la ruta del esfuerzo puramente matemático. Literalmente se corta y enmascara las propiedades de los campos estocásticos sin consultarle al solver implícito. 

**Anatomía de esta imposición**:
1. El integrador de `dustpy` recorre secuencias del `sim.gas.updateorder` (que se puede modificar usando `.insert(0, ...)`) .
2. Se fuerza matemáticamente un *Bump* en $\Sigma$ (e.g. `sim.gas.Sigma[:] = base * bump`).
3. El integrador corre, y efectivamente `P` (Presión del Disco) y $\rho$ (Densidad volumétrica) leen el salto brutal de $\Sigma$ exitosamente.
4. **Acoplamiento Perdido:** Para que haya un Gap u obstáculo para el polvo, se requiere que la **migración del arrastre de polvo reaccione a esto**. La velocidad radial del polvo depende paramétricamente de $\eta$ (el gradiente termodinámico de presión radial, definido numéricamente). 
   - Por razones de abstracción optimizada, en el marco estándar predeterminado de Tripod/Dustpy el actualizador analítico (no derivativo) de $\eta$ está divorciado del salto local e ignora la manipulación de base.
   - Si se prosigue por este método, es imperativo aplicar "Hacking en Cadena", donde además de hackear $\Sigma$, **hay que reemplazar `sim.gas.eta.updater.updater` por una función diferencial central puramente numérica** para forzar artificialmente que rastree las derivadas locales alteradas matemáticamente, evitando las suavizaciones analíticas subyacentes.

---

## Conclusión y Recomendación de Flujo:

Al analizar todo el marco expuesto, `dusty / tripodpy` brilla al extremo cuando se le permite solucionar sus ecuaciones de balance de manera puramente hidrodinámica (Modificación Subyacente).

1. **Para modelar vacíos hechos por planetas masivos (Gaps)**: El enfoque `dustpylib` (**Método 1**) de alterar $\alpha$ basándose en la prescripción analítica (Kanagawa/Duffell) o imponiendo un hueco propio en la viscosidad gaussiana es el estándar oro en la industria. Mantiene todos los balances termodinámicos.
2. **Para modelar sumideros / formación de cuerpos / zonas de sublimación**: Inyectar coeficientes de pérdida usando matrices negativas sobre **Método 2 (`dust.S.ext.updater`)**. 
3. **Para imponer *Muros / Bumps / Rampas* artificiales estrictamente definidas sin depender de decaimientos viscosos**: Habrá que sostener el **Método 3**. En cuyo caso hay que realizar un mapeo de derivadas manual (modificando numéricamente $\eta$ y si es conveniente, aislando derivativas en la torsión `v_rad`).

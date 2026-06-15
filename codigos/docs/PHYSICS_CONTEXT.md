# PHYSICS & PROJECT CONTEXT
> **Documento Maestro para IA y Desarrolladores**

Este archivo contiene el corpus teórico riguroso, el historial evolutivo, la arquitectura metodológica y las directrices de validación 1D del proyecto de acreción de pebbles y formación de "Waterworlds" en discos estructurados. Es de lectura obligatoria para cualquier IA al iniciar una sesión de trabajo en este repositorio.

## 1. Arquitectura Base: TripodPy y Simframe
El simulador subyacente es **TripodPy**, construido sobre la arquitectura de **Simframe**.
- **Grupos y Componentes:** Simframe organiza la simulación en `Groups` (`gas`, `dust`, `grid`). Cada grupo contiene `Fields` (ej. `gas.Sigma`).
- **Integración:** La evolución temporal se rige por una `IntegrationInstruction` que resuelve las ecuaciones de advección-difusión. La jerarquía estándar en TripodPy fuerza a que los parámetros de mezcla (`sim.dust.delta.rad`, `turb`, `vert`) hereden directamente el valor de `sim.gas.alpha`.
- *Nota arquitectónica:* Mantenemos consistencia importando constantes de `dustpy` y definimos grillas temporales logarítmicas en función de `c.year`. La salida genera archivos HDF5 por snapshot.

## 2. El Cambio de Paradigma: Pre y Post "Oka & Hartmann"
El proyecto sufrió una reestructuración masiva visible en el historial (`scratch/Pre oka and hartmann` vs `scratch/Post oka and hartmann`).
- **Pre-transición (Snowline estática):** Originalmente, la línea de nieve ($R_{snow}$) era una barrera fija. El flujo de pebbles cruzaba a una tasa constante determinada estrictamente por el drift radial.
- **Post-transición (Dinámica y Pebble Snow):** Se adoptó la física termo-viscosa donde la tasa de acreción estelar decae en el tiempo: $\dot{M}(t) = \dot{M}_{\text{1Myr}} \cdot (t/\text{1 Myr})^{-\eta}$. Siguiendo a **Oka et al. (2011)** y **Hartmann**, la posición de la snowline migra hacia el interior acoplada a este decaimiento (truncada entre ~2.73 AU y ~0.5 AU para estabilidad).
- **Justificación del Modelo Térmico:** Mantener esta interpolación secular analítica para $R_{snow}(t)$ elude la resolución de la ecuación de energía, lo que es óptimo y necesario para lograr cálculos estables y viables a $10$ Myr en régimen radiativo.

## 3. El Enfoque Teórico de Gijs Mulders: Arquitecturas y Volátiles
El entorno del modelo se fundamenta en los resultados del modelo PPOLs (McCloat, Mulders et al.).
- **Pebble Snow y Super-Tierras:** El flujo de guijarros de hielo cruzando la snowline empuja embriones hacia masas límite, derivando en Super-Tierras.
- **Waterworlds:** La posición dinámica de la snowline dicta la composición bulk. Formar un waterworld in situ depende del cruce entre la ventana de acreción y la migración secular.
- **Bimodalidad:** Dependiendo de la masa del disco, el sistema forma múltiples núcleos intra-snowline o bifurca a una arquitectura bimodal si los planetas gigantes externos alcanzan primero la masa de aislamiento.

## 4. Cinética del Polvo: Gradientes de Presión y Pressure Bumps
El polvo no interactúa gravitacionalmente con el disco, sino aerodinámicamente con el campo de presión del gas.
- **Gradiente de Presión ($\partial P / \partial r$) y Deriva:** El polvo sufre un headwind por orbitar en gas sub-kepleriano, migrando siempre hacia los máximos de presión.
- **Atrapamiento y Filtro Aerodinámico:**
  - *Trampas de Presión:* En $\partial P / \partial r = 0$, los sólidos convergen.
  - *Filtrado:* Un gap producido por $0.1 M_{Jup}$ con $\alpha=10^{-3}$ es un regulador ideal de flujo. Frena la deriva rápida atrapando pebbles de St marginal, pero no aísla totalmente el disco interior, permitiendo que el polvo fino ($\text{St} \ll \alpha$) fluya suavemente.

## 5. Física de Acreción de Pebbles (Ormel 2017 & Drążkowska et al. 2023)
Implementada en el módulo `PA3Py` (`PebbleAccretion3.py`).
- **Inicio de Pebble Accretion (PA Onset):** Ocurre cuando el decoplamiento gas-polvo ($St = \tau_f \Omega_K$) permite la captura gravitacional del embrión.
- **Masa de Aislamiento ($M_{iso}$):** El límite físico absoluto donde el embrión invierte el gradiente de presión ($\partial P / \partial r > 0$), anulando el flujo entrante de pebbles. Escala fuertemente con la razón de aspecto: $M_{iso} \propto (H/r)^3$.
- **Métrica de Éxito (Supervivencia):** Las simulaciones de acreción secular no buscan formar "Tierras" directas. El criterio de éxito físico es crear precursores masivos para la fase de impactos gigantes. **Todo núcleo que alcance o supere $\sim 0.1 M_\oplus$ (masa de Marte) es un éxito rotundo.**

## 6. Advertencias Críticas de Laboratorio vs Parametrización
- **La Disputa de $v_{frag}$ (Hielo vs Silicatos):** Tu modelo hace un salto de escalón en $v_{frag}$ (ej. 1000 cm/s para hielo a 100 cm/s para roca) para detonar el *Pebble Snow*. Sin embargo, **evidencia de laboratorio reciente dicta que el hielo fragmenta a velocidades casi idénticas a la roca** ($\sim 100-300$ cm/s).
- **Acción Obligatoria:** Es mandatorio correr casos de control donde $v_{frag}$ se mantenga **constante** (ej. 1 m/s) en todo el disco. Esto sirve para aislar el verdadero impacto aerodinámico del barrido de la snowline, separándolo del diferencial artificial de fragmentación térmica.

## 7. Pipeline y Extracción de Datos
- **Módulos Mixin (`pipeline_methods/`):** Estructura modular: `DiskSetupMixin`, `DiskChemistryMixin`, `SnowlinePhysicsMixin`, y `PressureBumpsMixin`.
- **Análisis Interactivo:** `analyzer_base.py` **siempre debe filtrar directorios por `alpha`** antes de inicializar HDF5/PA3Py para evitar corrupción de cachés (`.pkl`). Las gráficas en `plotters.py` siempre deben trazar la curva teórica de $M_{iso}(t)$.

## 8. Ejecución en Clúster (Creación de Jobs PBS)
El proyecto escala las simulaciones masivas (o la generación de visualizaciones como GIFs) mediante un sistema de colas en el **Clúster Geryon**.
- **Generación Dinámica de PBS:** Las instrucciones no se escriben a mano. Se usan scripts generadores (como `crear_pbs_gifs.py` o `lanzador_cluster_gifs.py`) que iteran mediante `glob` sobre los directorios de resultados (`data/runs/`) y construyen un archivo maestro (`.pbs`) concatenando comandos de python.
- **Entorno del Clúster:** Los encabezados PBS siempre deben incluir:
  - Carga de módulos de compilación: `module load devtoolset/devtoolset-11`.
  - Activación del entorno local de miniconda: `source /data4/maximiliano.valderrama/miniconda3/bin/activate tripod311`.
  - Asignación de recursos (ej. `#PBS -l select=1:ppn=8`, `#PBS -l walltime=24:00:00`).
- **Reemplazo de Rutas:** Cuando los generadores leen rutas en Windows (`\\`), el código siempre sanitiza y reemplaza por barras inclinadas (`/`) antes de inyectarlo al archivo PBS para evitar fallos de ruta en el clúster (Linux).

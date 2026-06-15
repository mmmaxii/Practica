# Conclusión y Toma de Decisiones: Generación de Subestructuras en TripodPy

## El Problema de los "Sinks" (Sumideros Manuales de Planetesimales)

A lo largo de nuestras experimentaciones con los métodos de la librería clásica de `dustpylib` (Drążkowska, Miller, Schoonenberg), logramos inyectar la física correctamente pero validamos que **el framework 1D de TripodPy (TwoPopPy) es inherentemente inestable para estas rutinas agresivas**, por dos motivos cruciales:

1. **Choque Dimensional Disfrazado:** El framework arrastra la aerodinámica y la densidad del plano medio ($\rho$, $St$, $a$) en 5 trazas (`shape = (150, 5)`), pero colapsa la masa integrada $\Sigma$ a 2 poblaciones. Modificar este comportamiento exige un *Fancy Indexing* adaptativo constante que rompe la simplicidad del código.
2. **Explosión Coagulativa Inocente ($10^9$ cm):** Al imponer una tasa veloz de captura de rocas (`zeta = 0.5`), el área elegida cae súbitamente en una baja densidad masiva (tendencia a cero). Dado que las ecuaciones de drift y velocidades por choques en un entorno 1D tienden a dividir por fracciones molares ligadas a la masa y a las tasas relativas al gas, la ecuación diverge. **El integrador entra en pánico y crea rocas de $10^9$ centímetros (10,000 km)**, destruyendo el realismo de la simulación. 

---

## La Solución "Estándar de Oro": Alteración Geodésica de la Viscosidad ($\alpha$)

Tras analizar la inestabilidad de remover masivamente polvo, se concluye oficialmente descartar el enfoque de "Sinks extractores" y centrar el proyecto íntegramente en la **Geodesia Fluidodinámica mediante perturbación de $\alpha$** (Kanagawa, Duffell o Bumps Gausianos Manuales).

**Ventajas Definitivas:**
1. **Paz Numérica:** En lugar de sustraer gramos matemáticamente, alteramos cómo fluye el gas. El solver implícito recalcula gentilmente presiones, temperatura y arrastre sin picos numéricos hostiles.
2. **Física Acoplada Hermosa (Pressure Traps):** Al hundir a $\alpha$ en ciertos puntos, el gradiente de presión reacciona cambiando el parámetro clave $\eta$. El *drift radial* empuja orgánicamente a los pebbles a surfear este gradiente, amontonándose limpiamente para crear anillos luminosos espectaculares sin divergencias en el $a_{max}$.
3. **El Truco del Salto en el Tiempo ($t=0$):** Al usar gaps analíticos estacionarios, podemos obviar esperar cientos de miles de órbitas e imponer la fisura instantáneamente adaptando el perfil en la inicialización:
   ```python
   sim.gas.Sigma[...] /= (sim.gas.alpha / alpha_base)
   sim.dust.Sigma[...] /= (sim.gas.alpha / alpha_base)[:, None]
   ```

---

## Directriz de Trabajo Futura

A partir de este diagnóstico, en adelante priorizaremos arquitecturas espaciales complejas (sistemas exoplanetarios tipo Trappist o Sistema Solar Joven) iterando la multiplicación de perfiles inversos sobre el vector basal de la turbulencia turbomagnética `sim.gas.alpha`.

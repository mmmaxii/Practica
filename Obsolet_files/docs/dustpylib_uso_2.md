# Diferencia Física en la Generación de Gaps y Trampas de Polvo

Al alterar el parámetro fundamental de viscosidad ($\alpha$) dentro de simulaciones hidro-termodinámicas como `dustpylib` / `tripodpy`, se pueden arrojar dos tipos de substructuras diametralmente opuestas según la prescripción utilizada. A continuación se resume la distinción física entre forzar una "Zona Muerta" a mano, versus delegar la tarea a una inyección de torque planetario de la librería estandarizada.

---

## 1. El Embudo Viscoso Manual (`alpha_dip_updater`) - "Trampas de Presión / Dead Zones"

Cuando se emplea una matemática manual (reducción gaussiana) de la viscosidad subyacente $\alpha$, se obliga violentamente al escalar de turbulencia (p.ej. de $10^{-3}$ a $10^{-4}$) en un tramo concéntrico preciso a paralizarse en lo profundo del disco. 

**¿Qué efecto físico induce sobre el framework?**
- **Efecto sobre el Gas:** Al decrecer de golpe la difusividad radial de la pista de gas (se crea la famosa "Dead Zone" o zona inactiva), el gas en lenta migración es incapaz de lidiar con un caudal superior de masa convergente desde afueras. Se incita el **"Efecto de Embotellamiento"**. Al no fluir rápido por la fricción baja, se acumula la densidad ($\Sigma_{gas}$) lo cual detona el alzamiento contundente de un verdadero Tensor Gravitacional-Presurizador (Pressure Bump natural).
- **Efecto sobre el Polvo:** Como las velocidades de fricción auto destructora molecular entre macro gránulos escalan linealmente acorde a factores de $\sqrt{\alpha}$, la violenta bajada de turbulencia los protege del destrozamiento por colisiones supersónicas. Adiciónale que la inalcanzable montaña de presión originada por el gas repele y detiene la mortal y constante Deriva Radial al pozo estelar. Obtiene como final una pista de polvo protegida e incontenida que crecerá hacia formas masivas de rocas y agregados milimétricos. Consumes tu simulación formando anillos constructivos inmensos y masivos.

## 2. El Integrador Analítico Planetario (`kanagawa2017`) - "Gaps Planetarios Destructivos"

Kanagawa simula en cambio un escenario donde se sitúa el poder de gravedad disolvente del torque expansivo inducido por un macro cuerpo orbital intruso (Plantea Gigante Júpiter-like). La función formula una ratio diferencial $f < 1$, que impone dictamen reverso.

**¿Qué efecto físico induce sobre el framework?**
- **Efecto sobre el Gas:** Establecer sobre el integrador iterativo la orden unánime de arrojar (`alpha0 / f`) provoca que tu base universal termodinámica de difusión explote agresivamente y simule el brutal y masivo torque del planeta repeliendo las trazas gaseosas perimetrales a sí mismo, tallando asintóticamente la caída local en vacío, es decir, el genuino "Planetary Gap". 
- **Efecto sobre el Polvo:** En las fronteras orbitarias del surco, la ausencia radical de contenciones gaseosas masivas lo sumerge a sufrir enormes índices de fricciosidad pura de la radiación basal de Alpha. Son inyectados al molino cósmico, devueltos de a poco y veloz a polvo residual o monomérico barriendose por los vientos o tragados al hoyo y estrellándose directo a su host star. Las únicas verdaderas repulsiones reconstructoras ocurren a miles de AUs, aglomerando gigantescas colmenas a modo de anillos bordeando los "labios o bordes agigantados de la represa de gas" o los picos que envuelven el agujero interno recién nacido. 

---

## Existen Otras Vías Oficiales de Generar Subestructuras

Dependiendo íntegramente del sustento físico original proyectado a modelizar: 

1. **Variabilidad Orgánica Planetaria Extensa (`duffell2020`):** 
   - Esta formulación evalúa variables parecidas pero bajo regímenes de números de Mach locales a las escalas de radio vector en interacción con órbitas. Operando casi de igual manera, este marco resulta ideal contra Kanagawa para pulir perfiles de modelación profundos contra someros, especialmente predilecto ante modelamiento de Planetas de perfiles Sub-Neptunianos.

2. **Sumideros Planetesimales A-Viscosos por Pebbles (`drazkowska2016`):**
   - Ruta revolucionaria en "Framework Sinks/Sumideros". Evades alterar turbulencias $\alpha$ ficticias y asumes total carencia en influencias por torque expansivo Planetario. Lo aplicas prescribiendo una coagulación espontánea tan letal por exceso denso y sobresaturación de las celdas macro-molares de Polvo/Gas que tu agregación cruza fronteras invisibles asumiéndose automáticamente **Rocas-Macro-Planetesimales**. 
   - Actúa borrando la métrica local y pasándola como coeficiente negativo observable desde el vector matriz base de Fuentes `sim.dust.S.ext.updater`. Como es un "efecto trampa", se chupa literalmente todo rastro visible de polvo causandos un gap visual fantasmal, impidiendo daños perjudiciales y maremotos a los dominios vecinos que no sienten nada.

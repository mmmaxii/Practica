# Pebble Accretion 3 (PA3Py)

Este módulo es una evolución algorítmica de `PebbleAccretion2.py` que integra estricta y directamente las formulaciones físicas obtenidas en:
- **Ormel (2017)**: *The Emerging Paradigm of Pebble Accretion*.
- **Drążkowska et al. (2023)**: *Planet Formation Theory in the Era of ALMA and Kepler: from Pebbles to Exoplanets*.

## Diferencias Clave con Versiones Anteriores

La arquitectura general del código (carga robusta y unificada de los HDF5 generados por `tripodpy` vía `from_datadir()`) se mantiene idéntica para asegurar total compatibilidad cruzada con el pipeline de evolución física, limitando exclusivamente los cambios a las fórmulas de régimen y cortes analíticos numéricamente evaluados.

### 1. Masa Crítica para iniciar la Acreción de Pebbles (*Onset Mass*)
Aparece una exclusión física para los cuerpos muy livianos. Para poder capturar el gas necesario de arrastre, un cuerpo debe tener una masa mínima. Implementamos:

$$ M_{\mathrm{PA\ onset}} = \mathrm{St}\, \eta^3 M_{\star} $$

*Implementación:* Si $M_{\rm pl} < M_{\rm PA\ onset}$, el modelo no asume tasa nula, sino que cruza al modo de **Acreción de Planetesimales (Safronov)** (Ec. 7.14 de Ormel 2017), asumiendo régimen gas-free balístico bajo densidad propia rocosa. Esto produce que haya una transición real entre el crecimiento primigenio de embriones pequeños al runaway estrepitoso provisto por pebbles.

### 2. Ecuaciones Dinámicas de Acreción
Se dividen orgánicamente los regímenes del acercamiento orbital según las propiedades fluidas analizadas por **Ormel (2017)**. 
- **Headwind (Bondi-like drift limit)** ($M \lesssim M_{\rm hw/sh}$):
$$ \dot{M}_{\mathrm{2D, hw}} = \sqrt{8 G M_{\mathrm{pl}} t_{\mathrm{stop}} v_{\mathrm{hw}}} \Sigma_{\mathrm{peb}} $$
- **Shear (Hill limit)** ($M \gtrsim M_{\rm hw/sh}$):
$$ \dot{M}_{\mathrm{2D, sh}} = 2 R_{\mathrm{Hill}}^2 \Omega_{\mathrm{K}} \mathrm{St}^{2/3} \Sigma_{\mathrm{peb}} $$
  
Con la masa de transición evaludaba matemáticamente con:
$$ M_{\mathrm{hw/sh}} = \frac{v_{\mathrm{hw}}^3}{8 G \Omega_{\mathrm{K}} t_{\mathrm{stop}}} $$

### 3. Transición Analítica a disco Extenso 3D
Ya no se evalúa `erf` numéricamente para forzar transiciones. Se utiliza la aproximación turbulenta del factor de decaimiento natural propuesta por **Dubrulle (1995) / Ormel 2017 (Eq 7.24)**:
$$ \dot{M} = \dot{M}_{\mathrm{2D}} \left( \frac{b_{\mathrm{col}}}{b_{\mathrm{col}} + h_{\mathrm{peb}} \sqrt{8/\pi}} \right) $$
Donde $h_{\rm peb} = \sqrt{\frac{\alpha_T}{\alpha_T + {\rm St}}} h_{\rm gas}$. Esta transición converge al límite 3D natural en discos muy agitados verticalmente.

### 4. Masa Reguladora de Aislamiento
Las barreras numéricas y el *overshoot* quedan protegidos por el rediseño más actualizado sobre retroalimentación entre aislamiento del gap y gas proporcionado por Bitsch, modificado bajo la última simplificación analítica de la comparación de **Drążkowska et al. 2023**:
$$ M_{\mathrm{iso, peb}} = 25 M_{\oplus} \left(\frac{H/r}{0.05}\right)^3 \left(\frac{M_{\star}}{M_{\odot}}\right) $$

### 5. Consumo y Drift
Los Pebbles viajan hacia la estrella limitados por:
$$ v_{\mathrm{r,solid}} \approx \frac{- 2\eta v_{\mathrm{K}} \mathrm{St}}{1 + \mathrm{St}^2} $$
Para mantener el módulo de acreción puramente parasitario y conservador, se extrae la traza dinámica de `tripodpy` asumiendo compatibilidad directa, validando los cruces del disco para compartir eficientemente los flujos con el acercamiento analógico del bucle espacial y temporal.

## Limitaciones y Notas Computacionales

- Para hacer cross-checking de perfiles puedes lanzar copias modulares intercambiando la clase originaria.
- El objeto instanciado expone un método llamado `.summary(results)` idéntico para que todo log del pipeline mantenga retrocompatibilidad.

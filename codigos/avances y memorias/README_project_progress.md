# Waterworld Formation via Pebble Accretion in Structured Disks
**Project Progress & Pipeline Documentation**

Este documento resume el avance total del proyecto de investigación sobre la formación de planetas tipo "Waterworld" mediante acreción de pebbles en discos protoplanetarios con subestructuras (gaps).

---

## 1. Visión General del Proyecto
El objetivo científico principal es investigar cómo las perturbaciones y subestructuras en un disco protoplanetario (gaps creados por planetas masivos o perturbaciones de viscosidad $\alpha$) afectan la eficiencia de la acreción de pebbles y la composición final de un embrión planetario situado en el disco interior. En particular, exploramos el mecanismo de **"pebble snow"**, donde el desplazamiento de la snowline libera un "burst" de pebbles helados retenidos en una trampa de presión.

### Mecanismos Físicos Implementados
1. **Dinámica del Disco (`tripodpy` / `dustpy`):** Evolución viscosa del gas y el polvo (drifting y coagulación/fragmentación), con seguimiento de especies químicas (Silicatos y H₂O).
2. **Modelo de Gaps:** Gaps implementados estáticamente al inicio de la simulación manipulando el perfil de viscosidad $\alpha(r)$ para evitar inestabilidades numéricas en el solver implícito.
3. **Pebble Accretion (PA3Py):** Integración de las masas planetarias utilizando el flujo de pebbles del disco. La composición se determina de manera dinámica al cruzar la snowline de H₂O.

---

## 2. Modelos y Subestructuras Estudiadas

Se desarrollaron tres arquitecturas principales de subestructuras, simuladas a $1 \times 10^6$ años (1 Myr) y $5 \times 10^6$ años (5 Myr):

- **Single Planet Gaps:** Un solo planeta masivo creando un gap. 
  - Masas: Super-Earth, Neptune, Saturn, Jupiter, Super-Jupiter.
  - Posiciones: 1.0, 2.0, 3.0, 5.0, 7.0 y 10.0 AU.
- **Multi-Planet Gaps:** Configuraciones de sistemas planetarios múltiples actuando como múltiples trampas de presión (ej. 4 Jupiters; Neptune + Saturn + Jupiter).
- **Sinusoidal Gaps:** Perturbaciones en el disco simulando oscilaciones de densidad a gran escala.
  - Configuraciones: Inner (1-20 AU) y Outer (5-50 AU) con 10 bumps equidistantes en escala logarítmica.
  - Amplitudes: 0.5, 1.0 y 2.0.

---

## 3. Arquitectura del Pipeline y Scripts

### 3.1 Pipeline de Simulaciones (`tripodpy`)
- `run_1myr.py`: Script central para ejecutar el grid completo de 37 simulaciones a 1 Myr. Incorpora refinamiento dinámico de la malla (resolución ajustada en los bordes de los gaps y la snowline) y guarda snapshots en HDF5.
- `run_5myr.py`: Extensión del pipeline para simulaciones a largo plazo (5 Myr). Actualmente en ejecución.
- `pipeline_methods/pressure_bumps.py`: Módulo "Mixin" que inyecta la física de gaps (Duffell, Sinusoidal, Gaussiano) modificando el parámetro $\alpha$ del disco.

### 3.2 Modelo de Acreción de Pebbles (`PA3Py`)
- `PebbleAccretion3.py`: Toma la densidad de gas, densidad de polvo, posiciones de snowlines y velocidades de deriva de los snapshots HDF5. Calcula la tasa de acreción $\dot{M}$ y divide la masa final entre silicatos, CO₂ y H₂O basándose en el estado termodinámico local (Ormel 2017, Drążkowska 2023).

### 3.3 Pipeline de Visualización y Análisis Científico
- `disk_visualizer/disk_viz.py`: Herramienta para visualizar en 2D la evolución morfológica del disco (gas y polvo). Genera animaciones estables fijando automáticamente las barras de escala logarítmica.
- `run_analysis_1myr.py`: Script automatizado para analizar lotes completos de simulaciones. Genera 4 figuras clave para investigación:
  1. **Heatmap f_H2O**: Mapa visual de la fracción final de H₂O en función de la masa y posición del gap.
  2. **Bar Chart Ordenado**: Clasifica las 37 runs de mayor a menor fracción de H₂O.
  3. **Evolución Temporal Agrupada**: Plots en escala log-log de la masa del embrión $M(t)$ divididos por categoría (Single, Multi, Sinusoidal).
  4. **Pebble Flux en la Snowline**: Diagnóstico directo del "Pebble Snow", mostrando $\dot{M}_{peb}$ cruzando la snowline.
- `run_analysis.bat`: Script en batch para Windows que automatiza la ejecución del análisis para embriones ubicados en 1.0, 2.0 y 3.0 AU secuencialmente y organiza los resultados en carpetas (`figs_1myr\rXau`).

---

## 4. Problemas Técnicos Resueltos Durante el Desarrollo

1. **Solver Hangs (Atascos del Solucionador):** Se corrigieron cuelgues del solver implícito al reducir la rigidez de las matrices. Los gaps se inyectaron de forma puramente estática (`setup_alpha_profile`) en $t=0$ para evitar saltos numéricos bruscos.
2. **Error "sim not defined":** Resolución de bugs de arquitectura en el pipeline modular (mixins) donde se perdía la referencia al objeto `Simulation` principal.
3. **Escalado Visual en Animaciones:** Corrección de parpadeos y desajustes de colorbars en `disk_viz.py` implementando persistencia de estados (`norms`) durante la renderización en MP4/GIF.
4. **Problemas con Entornos Virtuales (Unicode/Módulos):** Integración definitiva forzando UTF-8 en consola y garantizando que los scripts de análisis apunten al binario correcto de Python del entorno `env_tripod`.

---

## 5. Resultados Preliminares (Análisis de data_1myr)

Los primeros análisis para un embrión situado en 2.0 AU revelan hallazgos sumamente prometedores:
- **Formación Exitosa de Waterworlds:** 9 configuraciones lograron superar el umbral del 10% de H₂O en solo 1 Myr.
- **Single Gaps (El efecto de la posición):** 
  - Un Júpiter en 1.0 AU genera un asombroso **20.1% de H₂O**, actuando como regulador ideal que frena el polvo seco y deja pasar el "pebble snow" en el momento oportuno.
  - Planetas menos masivos (Super-Earth / Neptune) requieren estar en posiciones más externas (3 - 10 AU) para acumular ~13% de H₂O de manera eficiente.
- **Sinusoidal Gaps:** Una perturbación sinusoidal débil en el disco externo (`amp=0.5`) demuestra ser muy efectiva para atrapar hielo y liberarlo gradualmente (11.3% H₂O).
- **Múltiples Planetas:** Actúan como filtros en serie; múltiples gaps masivos bloquean casi totalmente el flujo hacia el disco interior, secando al embrión.

---

## 6. Siguientes Pasos
1. Concluir la ejecución del batch de **5 Myr** (`run_5myr.py`) que lleva en curso >30 horas.
2. Repetir el análisis `run_analysis_1myr.py` sobre los datos de 5 Myr para estudiar la convergencia a largo plazo de la composición del planeta.
3. Finalizar la generación de esquemas teóricos (visuales) usando IAs generativas de imagen basados en los metaprompts previamente diseñados.
4. Iniciar redacción de sección de resultados comparando la métrica de 1 Myr vs 5 Myr.

# Plan de Tesis: Modelado de Acreción de Pebbles en Discos Estructurados para estudiar la Formación de "Waterworlds"

En base a la propuesta de investigación y conversaciones con Gijs Mulders, el enfoque evolucionó hacia el uso de `tripodpy` como simulador de evolución del disco, aprovechando snowlines dinámicas y sub-estructuras de disco como los drivers naturales de trampas de presión. El modelo de acreción de pebbles (`PA3Py`) opera en post-proceso sobre los snapshots HDF5 del pipeline.

---

## Estado actual del pipeline (Abril 2026)

El pipeline `WaterworldPipeline` está completamente modularizado en el paquete `pipeline_methods/`:

| Módulo | Contenido |
|---|---|
| `disk_setup.py` | `setup_grid`, `setup_star`, `initialize_simulation`, `setup_star_evolution` |
| `disk_chemistry.py` | `add_volatile_components`, `add_ice_sigma_fields` |
| `snowline_physics.py` | `setup_physics` (v_frag(T)), `add_snowline_fields` (rsnow → HDF5) |
| `pressure_bumps.py` | `setup_gap_duffell`, `setup_gap_kanagawa`, `setup_alpha_sinusoidal`, `reset_gap` |

El módulo de acreción `PA3Py/PebbleAccretion3.py` lee los HDF5 y calcula la composición de los embriones a partir de `SigmaDust_X / Sigma_dust_total`.

---

## Fases de Implementación

### ✅ Fase 1: Química Multi-Especie y Snowlines Dinámicas
- Componentes volátiles ($H_2O$, $CO$, $CO_2$) inyectados como trazadores de tripodpy
- Perfil de $v_{\rm frag}(T)$ dinámico con escalones en cada snowline
- Contracción estelar pre-MS (Hayashi track) → migración de snowlines a lo largo de $10^7$ yr
- Fields `rsnow_X` y `SigmaDust_X` / `SigmaGas_X` guardados en cada snapshot HDF5

### ✅ Fase 2: Pipeline Maestro Modular
- Clase `WaterworldPipeline` con 4 mixins independientes
- Parámetros estelares ($M_\star, R_\star, T_\star$) configurables
- Output HDF5 gestionado por tripodpy con `diastole` updater para capturar estado post-solver

### ✅ Fase 3: Módulo de Acreción de Pebbles (PA3Py)
- Física de Ormel (2017) y Drążkowska et al. (2023): regímenes headwind/shear, masa de aislamiento
- Composición dinámica: $f_X = \Sigma_{\rm dust,X}(r,t) / \Sigma_{\rm dust,total}(r,t)$
- Lectura de `OmegaK` directamente del HDF5 (con fallback analítico)
- Trazado de $M_{\rm H_2O}$, $M_{\rm CO_2}$, $M_{\rm silicates}$ por embrión

### 🔄 Fase 4 (EN CURSO): Estudio del Efecto de Gaps en la Composición

**Objetivo científico:**
> ¿Cómo afectan los gaps y sub-estructuras del disco protoplanetario a la composición de los protoplanetas formados por acreción de pebbles?

**Modelo de gap:** Duffell (2020) — preferido sobre Kanagawa (2017) por su perfil empírico más suave y realista.

**Parámetros de simulación:**
- Integración: **10 Myr**, 50 snapshots log-uniformes
- Especie activa: **H₂O** (snowline principal)
- Grilla: 200 celdas, 1–200 AU, estrella solar (1 M☉)

**Grupos de simulaciones en `run_diagnostics.py`:**

#### Grupo 1 — Mismo tamaño, distintas posiciones
Un gap de 1 M♃ (Duffell) a 4 radios:

| Run | a_planet | M_planet |
|---|---|---|
| `duffell_1MJ_2au`  | 2 AU  | 1 M♃ |
| `duffell_1MJ_5au`  | 5 AU  | 1 M♃ |
| `duffell_1MJ_10au` | 10 AU | 1 M♃ |
| `duffell_1MJ_20au` | 20 AU | 1 M♃ |

#### Grupo 2 — Misma posición (5 AU), distintos tamaños
| Run | M_planet |
|---|---|
| `duffell_01MJ_5au` | 0.1 M♃ |
| `duffell_05MJ_5au` | 0.5 M♃ |
| `duffell_1MJ_5au`  | 1.0 M♃ |
| `duffell_2MJ_5au`  | 2.0 M♃ |

#### Grupo 3 — Múltiples gaps en distintas posiciones
| Run | Configuración |
|---|---|
| `duffell_multi_3a_5_7_10` | 1M♃@5AU + ½M♃@7AU + 1M♃@10AU |
| `duffell_multi_3b_sat_half_jup` | M♄@5AU + ½M♃@7AU + 1M♃@10AU |

#### Grupo 4 — Gaps sinusoidales (múltiples gaps equidistantes)
Perfil: $\alpha(r) = \alpha_0 \cdot [1 + A \cdot \sin^2(n \pi x(r))]$, con 5 gaps en 5–100 AU:

| Run | Amplitud A | Efecto |
|---|---|---|
| `sinusoidal_A1_suave`   | 1.0  | α_max ~ 2× baseline |
| `sinusoidal_A5_medio`   | 5.0  | α_max ~ 6× baseline |
| `sinusoidal_A10_fuerte` | 10.0 | α_max ~ 11× baseline |

### Fase 5 (PRÓXIMA): Comparativas Estelares y Demografía Planetaria
- Repetir Fase 4 para M-dwarfs (0.3–0.5 M☉) y estrellas A (1.5–2.0 M☉)
- Graficar distribuciones de masa vs. fracción de agua por tipo estelar
- Identificar el espacio de parámetros de disco que produce waterworlds en la zona habitable conservadora

---

## Pregunta científica central

> Para una estrella solar con distintas arquitecturas de sub-estructura en el disco, ¿en qué rango de masas de disco y posiciones de gap se maximiza la entrega de hielo ($H_2O$) a embriones en la zona habitable (~1 AU), produciendo planetas con fracción de agua Earth-like (~0.01–0.1%)?

---

## Archivos clave

| Archivo | Función |
|---|---|
| `codigos/pipeline_snowlines.py` | Pipeline principal (clase WaterworldPipeline) |
| `codigos/pipeline_methods/` | Módulos mixin (setup, química, snowlines, gaps) |
| `codigos/run_diagnostics.py` | Script multi-run con todas las configuraciones |
| `codigos/plot_diagnostics.py` | SnowlineDiagnostics (6 figuras por run) |
| `codigos/PA3Py/PebbleAccretion3.py` | Módulo de acreción (Ormel 2017 + composición) |
| `codigos/PA3Py/plot_pa3.py` | Visualización de resultados de acreción |

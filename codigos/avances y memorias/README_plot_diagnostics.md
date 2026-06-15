# `plot_diagnostics.py` — Guía de uso

Herramienta de diagnóstico visual para simulaciones de **tripodpy** con snowlines y acumulación de pebbles.  
Lee directamente los archivos HDF5 generados por el pipeline y produce los siguientes gráficos:

1. Hovmöller `(eps / Sigma_dust / Sigma_gas)` — diagrama espacio-tiempo
2. Hovmöller R(t) estilo artículo — multipanel con radio en Y, tiempo en X
3. Distribución de tamaño de grano — Σ(a, r) en una snapshot
4. Pebble flux Ṁ_pebble(r, t) — flujo másico espacio-temporal
5. Perfiles de disco (η, St, a_max, Σ) — 4 paneles en un snapshot
6. Densidades superficiales de gas por componente
7. Densidades superficiales de polvo por componente

---

## Uso

```bash
# Desde terminal (con el env activo: env_tripod)
python plot_diagnostics.py <datadir> [savedir]

# Ejemplo
python plot_diagnostics.py output_test_pipeline diagnostics_output/
```

```python
# Como módulo en un notebook
from plot_diagnostics import SnowlineDiagnostics

d = SnowlineDiagnostics(
    "output_test_pipeline",
    savedir="figs",
    r_trim=0.93,       # Conserva el 93% interior del grid radial
    t_min_yr=100.0     # Omite snapshots antes de t=100 años
)
d.plot_hovmoller(quantity="eps")
d.plot_hovmoller(quantity="Sigma_dust")
d.plot_hovmoller_rt(quantities=["Sigma_gas", "T", "a_max"], t_unit="kyr")
d.plot_size_distribution(it=-1)
d.plot_pebble_flux()
d.plot_profiles(it=-1)
d.plot_gas_components(it=-1)
d.plot_dust_components(it=-1)
plt.show()
```

---

## Parámetros de configuración principal (`SnowlineDiagnostics`)

| Parámetro | Default | Descripción |
|---|---|---|
| `datadir` | — | Directorio con los archivos `data0000.hdf5`, `data0001.hdf5`... |
| `savedir` | `None` | Directorio donde guardar PDFs. Si es `None`, no guarda. |
| `r_trim` | `0.93` | Fracción del grid radial a conservar. Recorta el borde exterior para evitar artefactos del `SigmaFloor` (explosiones numéricas). |
| `t_min_yr` | `100.0` | Snapshot temporal mínimo en años. Omite los primeros pasos donde el disco no tiene polvo físicamente relevante. |

---

## Gráficos generados

### 1. `plot_hovmoller(quantity)` — Diagrama espacio-tiempo (Hovmöller)
- **Ejes**: Tiempo [años] (X, log) vs Radio [AU] (Y, log)
- **Color**: `eps` = relación polvo/gas ε | `Sigma_dust` | `Sigma_gas`
- **Colormap**: `magma`
- Superpone curvas `r_snow(t)` para las snowlines de H₂O, CO₂ y CO

| Parámetro | Default | Descripción |
|---|---|---|
| `quantity` | `"eps"` | Cantidad a graficar: `"eps"`, `"Sigma_dust"`, `"Sigma_gas"` |
| `cmap` | `"magma"` | Colormap de matplotlib |
| `vmin`, `vmax` | `None` | Límites del colorbar (log₁₀). Si `None`, auto-escala por percentiles. |

---

### 1b. `plot_hovmoller_rt(quantities, ...)` — Hovmöller multi-panel estilo artículo
- **Ejes**: Tiempo (X) vs Radio [AU] (Y, log) — formato estándar de papers (Drążkowska et al., Booth & Ilee, etc.)
- Genera N paneles apilados verticalmente con eje X compartido
- Permite ver simultáneamente la evolución de Σ_gas, T, a_max, ε, y St
- Snowlines superpuestas como curvas `r_snow(t)` en cada panel

| Parámetro | Default | Descripción |
|---|---|---|
| `quantities` | `["Sigma_gas", "T", "a_max"]` | Lista de cantidades a mostrar. Opciones: `"Sigma_gas"`, `"Sigma_dust"`, `"T"`, `"a_max"`, `"eps"`, `"St"` |
| `cmaps` | `None` | Dict con colormaps por cantidad. Ej: `{"T": "inferno"}` |
| `t_unit` | `"kyr"` | Unidades del eje temporal: `"yr"`, `"kyr"`, `"Myr"` |
| `t_log` | `False` | Si `True`, eje X en escala logarítmica |
| `figsize` | auto | Tupla. Default: `(11, 3.5 × n_paneles)` |
| `percentile_range` | `(2, 98)` | Percentiles para la escala de color de cada panel |

---

### 2. `plot_size_distribution(it)` — Distribución de tamaño de grano
- **Ejes**: Radio [AU] (X, log) vs Tamaño `a` [cm] (Y, log)
- **Color**: Densidad superficial Σ(a, r) [g/cm²] (log)
- TriPoD usa 2 bins (`dust.Sigma`); reconstruye el mapa distribuyendo cada bin como gaussiana en log(a) centrada en `s.min` (bin 0) y `s.max` (bin 1)
- Superpone `a_max(r)` y `a_min(r)` como líneas blancas

| Parámetro | Default | Descripción |
|---|---|---|
| `it` | `-1` | Índice temporal dentro de los snapshots válidos (-1 = último) |
| `cmap` | `"viridis"` | Colormap |
| `vmin`, `vmax` | `None` | Límites del colorbar en log₁₀ |
| `sigma_spread` | `0.5` | Ancho en décadas de log(a) de la gaussiana por bin |

---

### 3. `plot_pebble_flux()` — Flujo másico de pebbles
- **Ejes**: Tiempo [años] (X, log) vs Radio [AU] (Y, log)
- **Color**: |Ṁ_pebble| en M⊕/yr (log₁₀), colormap `RdYlBu_r`
- Calcula Ṁ = −2π r Σ_dust v_r usando velocidad ponderada por masa:
  - bin 0 (pequeño): `v.rad[..., 2]` (velocidad en `a1`)
  - bin 1 (grande): `v.rad[..., 4]` (velocidad en `amax`)
- `dust.v.rad` tiene shape `(Nt, Nr, 5)` con velocidades en `[a0, fudge·a1, a1, 0.5·amax, amax]`
- Superpone curvas `r_snow(t)` para las snowlines

| Parámetro | Default | Descripción |
|---|---|---|
| `cmap` | `"RdYlBu_r"` | Colormap |
| `vmin`, `vmax` | `None` | Límites del colorbar en log₁₀ |

---

### 4. `plot_profiles(it)` — Perfiles físicos del disco
4 subpaneles en una snapshot:
- **η**: parámetro de gradiente de presión (panel superior izquierdo)
- **St**: número de Stokes del bin más grande (panel superior derecho)
- **a_max / a_min**: tamaños de grano límite (panel inferior izquierdo)
- **Σ_gas + Σ_dust**: densidades superficiales totales (panel inferior derecho)

Snowlines marcadas como líneas verticales en la posición exacta del snapshot indicado.

| Parámetro | Default | Descripción |
|---|---|---|
| `it` | `-1` | Índice del snapshot (-1 = último) |

---

### 5. `plot_gas_components(it)` — Densidades de gas por componente
- Σ_gas de cada componente volátil activo (H₂O, CO₂, CO, Default...) leído desde el dump file
- Gas total (desde HDF5) en línea negra discontinua
- Leyenda exterior para mayor claridad
- Snowlines como líneas verticales en la posición correcta para ese snapshot

| Parámetro | Default | Descripción |
|---|---|---|
| `it` | `-1` | Índice del snapshot para el título (-1 = último) |
| `ylim` | `(1e-5, 1e3)` | Límites del eje Y |
| `legend_outside` | `True` | Si `True`, coloca la leyenda a la derecha del gráfico |
| `dumpfile` | `None` | Ruta al `frame.dmp`. Si `None`, busca en `datadir/frame.dmp` |

---

### 6. `plot_dust_components(it)` — Densidades de polvo por componente
- Σ_dust de cada componente sólido activo (H₂O, CO₂, CO, silicatos...) leído desde el dump file
- Silicatos graficados en negro al final como referencia
- Snowlines como líneas verticales en la posición correcta para ese snapshot

| Parámetro | Default | Descripción |
|---|---|---|
| `it` | `-1` | Índice del snapshot para el título (-1 = último) |
| `ylim` | `(1e-8, 1e3)` | Límites del eje Y |
| `legend_outside` | `True` | Si `True`, coloca la leyenda a la derecha del gráfico |
| `dumpfile` | `None` | Ruta al `frame.dmp`. Si `None`, busca en `datadir/frame.dmp` |

---

## ⚠️ Requisito para los plots de componentes

Las densidades superficiales por componente **no se guardan en los archivos HDF5 por defecto**. Los métodos `plot_gas_components` y `plot_dust_components` leen los datos desde el **dump file `frame.dmp`**, que sí contiene el estado completo de la simulación (serializado con `dill`).

El `frame.dmp` se genera automáticamente cuando la simulación hace pause o checkpoint. Si no existe, el script imprime instrucciones para activar `save = True` manualmente:

```python
# Alternativa si no hay frame.dmp: activar guardado en los HDF5
sim.components.H2O.gas.Sigma.save  = True
sim.components.H2O.dust.Sigma.save = True
sim.components.CO2.gas.Sigma.save  = True
# etc.
```

> Si no hay `frame.dmp` ni guardado activado, los gráficos de componentes mostrarán solo el gas total.

---

## Snowlines dinámicas

Las snowlines no se marcan siempre en la posición del primer snapshot (t=0). El sistema usa la posición correcta para cada tiempo:

| Tipo de plot | Método | Descripción |
|---|---|---|
| Hovmöller, Pebble flux (2D) | `_add_snowline_curves()` | Curva `r_snow(t)` completa superpuesta al mapa de color |
| Perfiles, Distribución de tamaño, Componentes (snapshot) | `_add_snowline_vlines(it_abs=...)` | Línea vertical en la posición exacta del snowline en ese tiempo |

### Cómo se leen las snowlines

**Prioridad 1 — Desde el Field guardado en HDF5** (`add_snowline_fields()` en el pipeline):

```python
sim.grid.addfield("rsnow_H2O", 0., description="Snowline H2O [cm]")

def rsnow_updater(sim):
    isnow = np.argmax(sim.gas.T < T_sub)
    return float(sim.grid.ri[isnow])

sim.grid.rsnow_H2O.updater.updater = rsnow_updater
```

Como son `Field`s de simframe con `save=True` (default), se guardan automáticamente en cada HDF5. Después de `read.all()` tienen shape `(Nt,)`.

El sistema incluye una validación de sanidad: si ≥ 80% de los valores están pegados al borde interior (bug de `T_bind` en vez de `T_sub`), descarta el Field y cae al fallback.

**Fallback — Cálculo desde `gas.T`**:

Si los HDF5 fueron generados sin `add_snowline_fields()`, el sistema calcula las snowlines snapshot a snapshot desde el perfil de temperatura, sin necesidad de intervención manual.

---

## Diagnóstico

Si los gráficos de componentes aparecen vacíos, usar el método `debug_components()` para inspeccionar qué está disponible en los datos leídos del HDF5:

```python
from plot_diagnostics import SnowlineDiagnostics
d = SnowlineDiagnostics("output_test_pipeline")
d.debug_components()
```

Imprime la estructura completa de cada componente (atributos, shapes, etc.) tal como quedaron guardados en el HDF5.

---

## Detalles técnicos

- **Unidades internas**: totalmente CGS (como las almacena tripodpy)
- **Ejes de gráficos**: AU y años (o kyr/Myr) para legibilidad física
- **Formato de salida**: únicamente PDF (alta calidad vectorial)
- **Snowline detection**: busca donde `T < T_sub` para H₂O (150 K), CO₂ (70 K), CO (25 K)
- **Lector de datos**: usa `dustpy.hdf5writer().read.all()` para leer todos los snapshots
- **Componentes**: leídos desde `frame.dmp` con `dill` (pickle extendido)

### Colores predefinidos por componente

| Componente | Color hex |
|---|---|
| H₂O | `#4FC3F7` |
| CO₂ | `#80CBC4` |
| CO | `#FFCC80` |
| silicates | `#EF9A9A` |
| Default | `#CE93D8` |

---

## Estructura del datadir esperada

```
output_test_pipeline/
    data0000.hdf5    ← snapshot t=0
    data0001.hdf5    ← snapshot t=...
    ...
    data0019.hdf5    ← snapshot final
    frame.dmp        ← dump file completo (necesario para plot_*_components)
```

### Campos guardados en los HDF5 por defecto

| Campo | Shape | Guardado por defecto |
|---|---|---|
| `gas/Sigma` | (Nr,) | ✅ Sí |
| `gas/T` | (Nr,) | ✅ Sí |
| `gas/eta` | (Nr,) | ✅ Sí |
| `dust/Sigma` | (Nr, 2) | ✅ Sí |
| `dust/v/rad` | (Nr, 5) | ✅ Sí |
| `dust/St` | (Nr, 5) | ✅ Sí |
| `dust/s/max`, `dust/s/min` | (Nr,) | ✅ Sí |
| `grid/rsnow_H2O`, `grid/rsnow_CO2`, `grid/rsnow_CO` | (Nt,) | ✅ Si se llama `add_snowline_fields()` |
| `components/H2O/gas/Sigma` | (Nr,) | ❌ Solo en `frame.dmp` |
| `components/H2O/dust/Sigma` | (Nr, 2) | ❌ Solo en `frame.dmp` |

---

## Notas de física

- **Recorte radial `r_trim`**: El borde exterior del disco cae al `SigmaFloor` causando artefactos numéricos en η, St, y Σ. Se recomienda `r_trim=0.90–0.95`.
- **Recorte temporal `t_min_yr`**: Los primeros snapshots (< 100 yr) no tienen polvo apreciable. El flujo másico y el diagrama Hovmöller son físicamente irrelevantes en ese rango.
- **Distribución de tamaño**: TriPoD (Two-Population Dust) usa exactamente 2 bins de Sigma, a diferencia de DustPy que usa muchos bins. El mapa 2D en log(a) es una reconstrucción aproximada mediante gaussianas de ancho configurable (`sigma_spread`).
- **Velocidades de fragmentación**: el pipeline implementa `v_frag` variable con la temperatura: 100 cm/s (silicatos secos), 1000 cm/s (hielo de H₂O), 500 cm/s (hielo de CO). Esto genera el "traffic jam" en las snowlines.
- **`plot_hovmoller_rt`**: formato estándar de artículos (Drążkowska et al., Booth & Ilee, etc.) con radio en eje Y. Permite ver simultáneamente cómo evolucionan estructuras radiales (snowlines, bumps, traffic jams) a lo largo del tiempo.

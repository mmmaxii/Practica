# make_gifs_5myr.py — Generación batch de GIFs para data_5myr

Genera automáticamente animaciones GIF del disco protoplanetario para **todos**
los runs de `data_5myr/`, llamando a `disk_visualizer/disk_viz.py` por cada uno.
Cada GIF se guarda dentro del propio directorio del run, en una subcarpeta `gif/`.

---

## Estructura de salida

```
data_5myr/
├── single_jup_3.0au/
│   ├── data0000.hdf5
│   ├── ...
│   ├── frame.dmp
│   └── gif/
│       └── disco.gif        ← generado aquí
├── single_sat_5.0au/
│   └── gif/
│       └── disco.gif
└── ...
```

---

## Uso básico

```powershell
# Con los defaults (every=2, fps=5, rmax=30 AU)
py make_gifs_5myr.py

# Personalizar parámetros
py make_gifs_5myr.py --every 3 --fps 8 --rmax 25

# Omitir GIFs que ya existen (ideal para reanudar una ejecución interrumpida)
py make_gifs_5myr.py --skip-existing

# Procesar solo algunos runs específicos
py make_gifs_5myr.py --runs single_jup_3.0au single_sat_5.0au multi_jup5_hjup7_jup10

# Usar una carpeta de datos diferente
py make_gifs_5myr.py --datadir data_1myr
```

---

## Referencia de argumentos

| Argumento | Tipo | Default | Descripción |
|-----------|------|---------|-------------|
| `--datadir` | `str` | `data_5myr` | Carpeta raíz con los runs |
| `--every` | `int` | `2` | Usar 1 de cada N snapshots |
| `--fps` | `int` | `5` | Frames por segundo del GIF |
| `--rmax` | `float` | `30` | Radio máximo del plot \[AU\] |
| `--dpi` | `int` | `100` | DPI de los frames renderizados |
| `--gif-name` | `str` | `disco.gif` | Nombre del archivo GIF dentro de `gif/` |
| `--runs` | `str...` | todos | Lista de runs a procesar (nombres de carpeta) |
| `--skip-existing` | flag | `False` | Omitir runs cuyo GIF ya existe |

---

## Comportamiento

- **Detección automática de runs**: encuentra todos los subdirectorios de
  `--datadir` que contengan al menos un archivo `data*.hdf5`.
- **Runs sin HDF5**: se listan como `[!] omitido` pero no interrumpen el batch.
- **Errores por run**: si `disk_viz.py` falla en un run, el script continúa
  con los demás y reporta los fallos al final.
- **Resumen final**: imprime cuántos GIFs se generaron correctamente y cuáles fallaron.

---

## Ejemplo de salida en consola

```
Buscando runs en: data_5myr/
  [!] _tmp: sin HDF5 — omitido

======================================================================
  Runs a procesar : 35
  every=2  fps=5  rmax=30.0 AU  dpi=100
  GIF → <run>/gif/disco.gif
  skip-existing   : False
======================================================================

--- [1/35] single_jup_1.0au ---
======================================================================
[GIF] single_jup_1.0au
  → Salida: data_5myr/single_jup_1.0au/gif/disco.gif
======================================================================
Generando GIF con 30 frames → data_5myr/single_jup_1.0au/gif/disco.gif
  ✓ Completado en 1.4 min  |  3.2 MB

...

======================================================================
RESUMEN: 35/35 GIFs generados correctamente
======================================================================
```

---

## Dependencias

Requiere las mismas dependencias que `disk_viz.py`:

```powershell
pip install h5py numpy matplotlib Pillow tqdm
```

---

## Ver también

- [`disk_visualizer/README.md`](disk_visualizer/README.md) — documentación completa de `disk_viz.py`
  (modos interactivo, snapshot único, frames PNG, referencia de argumentos)

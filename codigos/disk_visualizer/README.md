# disk_viz.py — Visualizador 2D de Disco Protoplanetario

Genera animaciones y capturas 2D de un disco protoplanetario a partir de los
snapshots HDF5 producidos por **tripodpy**. Los perfiles radiales 1D se
"revuelven" azimutalmente para producir una imagen tipo cara-frontal del disco.

---

## Paneles visualizados

```
┌──────────────────────┬──────────────────────┐
│   Σ_gas  [g cm⁻²]   │   Σ_dust [g cm⁻²]   │
└──────────────────────┴──────────────────────┘
```

- La **snowline H₂O** se dibuja como un contorno blanco punteado con su radio en AU.
- La **estrella central** aparece como un punto brillante en el origen.
- Escalas de color logarítmicas con rango fijo en el primer frame (colorbars estáticas).

---

## Instalación de dependencias

```powershell
pip install h5py numpy matplotlib Pillow tqdm
```

---

## Modos de uso

### 1. GIF animado *(modo principal)*

```powershell
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au --gif disco.gif
```

Con parámetros personalizados:

```powershell
py disk_viz.py `
    --datadir ../data_5myr/single_jup_3.0au `
    --gif     disco_jup3.gif `
    --every   2 `
    --fps     5 `
    --rmax    30 `
    --dpi     100
```

---

### 2. Frames PNG individuales

Guarda cada snapshot como `frame_XXXX.png` en la carpeta de salida:

```powershell
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au --output frames/
```

---

### 3. Snapshot único

Renderiza y guarda un único frame por índice (0-based):

```powershell
# Primer snapshot
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au --snapshot 0

# Último snapshot (índice negativo no disponible por CLI — usar índice real)
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au --snapshot 59
```

---

### 4. Modo interactivo *(requiere display)*

Navega por los snapshots con las teclas de flecha:

```powershell
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au --interactive
```

| Tecla | Acción |
|-------|--------|
| `→` / `d` / `Espacio` | Siguiente snapshot |
| `←` / `a` | Snapshot anterior |
| `q` | Cerrar |

---

### 5. Preview rápida *(sin argumentos extra)*

Si no se pasa `--gif`, `--output`, `--snapshot` ni `--interactive`, guarda
`disco_inicio.png` y `disco_final.png` en el directorio actual:

```powershell
py disk_viz.py --datadir ../data_5myr/single_jup_3.0au
```

---

## Referencia de argumentos

| Argumento | Tipo | Default | Descripción |
|-----------|------|---------|-------------|
| `--datadir` / `-d` | `str` | *requerido* | Directorio con los HDF5 del run |
| `--gif` / `-g` | `str` | `None` | Ruta del GIF de salida |
| `--output` / `-o` | `str` | `None` | Carpeta para guardar frames PNG |
| `--interactive` / `-i` | flag | `False` | Modo interactivo con teclado |
| `--snapshot` / `-s` | `int` | `None` | Índice de snapshot único a guardar |
| `--rmax` | `float` | `40` | Radio máximo del plot \[AU\] |
| `--every` | `int` | `1` | Usar 1 de cada N snapshots |
| `--fps` | `int` | `6` | Frames por segundo del GIF |
| `--dpi` | `int` | `120` | DPI de los frames renderizados |

---

## Campos HDF5 leídos

| Campo HDF5 | Descripción |
|------------|-------------|
| `grid/r` | Radio radial \[cm\] |
| `t` | Tiempo de simulación \[s\] |
| `gas/Sigma` | Densidad superficial de gas \[g cm⁻²\] |
| `gas/T` | Temperatura del gas \[K\] |
| `gas/alpha` | Parámetro de viscosidad α |
| `dust/Sigma` | Densidad superficial de polvo (sumada sobre bins) |
| `grid/SigmaDust_H2O` | Σ de polvo H₂O (campo personalizado de tripodpy) |
| `grid/rsnow_H2O` | Radio de snowline H₂O \[cm\] |

---

## Batch para todos los runs de data_5myr

Usar el script `make_gifs_5myr.py` en la carpeta padre:

```powershell
cd ..
py make_gifs_5myr.py --every 2 --fps 5 --rmax 30
```

Cada GIF se guarda en `data_5myr/<run>/gif/disco.gif`.
Ver `README_make_gifs.md` para más detalles.

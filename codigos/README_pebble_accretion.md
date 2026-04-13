# README — Módulo de Acreción de Pebbles (`PebbleAccretion2.py`)

> **Estado:** Funcional y validado con datos reales de `tripodpy` (pipeline_icefrac, 2026-04-06)  
> **Autores:** Maxi + Antigravity (desarrollado iterativamente desde `PebbleAccretion.py` → `v2` → corregido)

---

## 1. Resultados de la prueba de integración

El módulo fue probado con los datos de `data_post_pipeline/pipeline_icefrac/`:
- **Grid:** 1 – 1000 AU, Nr = 300, Nt = 30 snapshots
- **Tiempo simulado:** t = 10³ – 10⁵ años
- **Estrella:** contrae de 3 R☉ → 1 R☉ en 10⁴ años (snowlines migran hacia adentro)

### Resultados por embrión (M₀ = 0.012 M⊕ ≈ 1 M_luna)

| r [AU] | M_final [M⊕] | H₂O% | CO₂% | Sil% | Tipo |
|-------:|-------------:|-----:|-----:|-----:|------|
| 1.0 | 9.46 | 0.0 | 0.0 | 99.9 | Rocoso |
| 2.0 | 13.2 | 0.0 | 0.0 | 99.9 | Rocoso |
| 3.0 | 5.90 | 0.0 | 0.0 | 99.8 | Rocoso |
| **4.0** | **6.77** | **7.4** | 0.0 | 92.4 | Rocoso |
| **5.0** | **6.00** | **7.4** | 0.0 | 92.4 | Rocoso |
| **8.0** | **3.23** | **7.4** | 0.0 | 92.2 | Rocoso |
| 15.0 | 0.052 | 5.7 | 0.15 | 70.7 | Rocoso |
| 20.0 | 0.018 | 2.4 | 0.6 | 29.7 | Rocoso |
| 30.0 | 0.013 | 0.5 | 0.1 | 6.7 | Rocoso |

**Snowlines en t_final:**
- rsnow_H₂O = 3.24 AU
- rsnow_CO₂ = 15.14 AU
- rsnow_CO = 117.49 AU

### Interpretación física

- **≤3 AU** → 100% silicatos: embriones dentro del snowline H₂O durante toda la simulación → sin hielo ✅
- **4–8 AU** → ~7.4% H₂O: fuera del snowline (rsnow_H₂O inicialmente en ~30 AU con R★=3R☉, luego migra a 3.24 AU). Al principio sí había hielo disponible. ✅
- **>15 AU** → H₂O cae porque a t=0 el snowline todavía estaba ahí (disco caliente), y el dust drift dejó poco material disponible a radios grandes.
- **Ninguno supera 10% H₂O** porque el tiempo de simulación (10⁵ yr) es corto: la estrella contrae en 10⁴ yr pero los embriones acumulan composición en una ventana temporal estrecha justo después del cruce del snowline.

> Para obtener Waterworlds (f_H₂O > 10%) se necesita: simulaciones hasta ~1 Myr, embriones plantados en 4–10 AU, y/o masa inicial menor (semillas más primitivas).

---

## 2. Arquitectura del sistema

```
pipeline_snowlines.py
  ├── add_volatile_components()        ← inyecta H2O, CO2, CO (y más) desde chem.txt
  │     active_species = ["H2O","CO2","CO"]   ← configurable antes de correr
  ├── setup_physics()                  ← v_frag(T) dinámico según active_species
  ├── add_snowline_fields()            ← guarda grid/rsnow_X en cada HDF5
  ├── add_ice_sigma_fields()           ← guarda grid/SigmaDust_X y grid/SigmaGas_X
  │     Updater real: comp.dust.Sigma.sum(-1)  → shape (Nr,) por snapshot
  └── run_integration()               → data_post_pipeline/*/  (Nt HDF5)

PebbleAccretion2.py
  └── PebbleAccretionModule
        ├── from_datadir(datadir)     ← apila Nt HDF5 en arrays (Nt, Nr, ...)
        │     Auto-detecta grid/SigmaDust_X  →  volatile species
        │     Silicatos = dust_total − Σ(SigmaDust_X)  [físicamente correcto]
        │     Prioridad: SigmaDust_X > SigmaIce_X > components/*/dust/Sigma
        ├── run_growth(embryos_AU)    ← integra dM/dt con clamp M_iso por paso
        │     Composición: frac_X = SigmaDust_X / Σ_Y SigmaDust_Y  [real]
        └── summary(results)          ← tabla resumen

test_pebble_accretion.py             ← script de prueba, genera 3 PDFs en figs_pebbles/
diag_components.py                   ← script de diagnóstico numérico
```

---

## 3. Física implementada y referencias

### 3.1 Velocidad relativa pebble–embrión

$$\Delta v = \eta \, v_K, \qquad v_K = \Omega_K \cdot r$$

> Lambrechts & Johansen (2012), A&A 544, A32, **Eq. 16**

### 3.2 Radio de captura efectivo (settling regime)

$$R_{\rm acc} = \min\left(R_H,\ b_{\rm Bondi},\ b_{\rm Hill}\right)$$

$$R_H = r\left(\frac{M}{3M_*}\right)^{1/3}, \quad b_{\rm Bondi} = \sqrt{St}\frac{GM}{\Delta v^2}, \quad b_{\rm Hill} = \left(\frac{St}{0.1}\right)^{1/3}R_H$$

> Ormel & Klahr (2010), A&A 520, A43, Eq. A9 — Lambrechts & Johansen (2012), Eq. A10

### 3.3 Altura del disco de pebbles

$$H_{\rm peb} = H_{\rm gas}\sqrt{\frac{\alpha}{\alpha + St}}$$

> Youdin & Lithwick (2007), Icarus 192, 588, Eq. 12

### 3.4 Tasas de acreción y transición 2D/3Dxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx·

$$\dot{M}_{2D} = 2\,R_{\rm acc}\,\Sigma_{\rm peb}\,\Delta v, \qquad \dot{M}_{3D} = \pi\,R_{\rm acc}^2\,\rho_{\rm peb}\,\Delta v$$

$$\dot{M}_{\rm core} = \mathrm{erf}(z)\,\dot{M}_{2D} + [1-\mathrm{erf}(z)]\,\dot{M}_{3D}, \quad z = \frac{R_{\rm acc}}{\sqrt{2}\,H_{\rm peb}}$$

> Liu & Ormel (2018), A&A 615, A178, **Appendix A.1**

### 3.5 Masa de aislamiento

$$M_{\rm iso} = 25\,M_\oplus\left(\frac{h}{0.05}\right)^3 f_\alpha\,f_P$$

$$f_\alpha = 0.34\left[\log_{10}\!\left(\frac{\alpha}{10^{-3}}\right)+3\right]+0.66, \quad f_P = 1 - \frac{\partial\ln P/\partial\ln r + 2.5}{6}$$

$$\frac{\partial\ln P}{\partial\ln r} = \frac{-2\eta}{(H/r)^2}$$

> Bitsch et al. (2018), A&A 612, A30, Eq. 5 — Lambrechts et al. (2014), A&A 572, A35 (fallback)

### 3.6 Composición química — campos reales por componente

**Campo guardado en cada snapshot HDF5 (método actual):**

```
grid/SigmaDust_{sp}(r, t)  ←  comp.dust.Sigma.sum(-1)   [g/cm²]  densidad real de polvo sólido
grid/SigmaGas_{sp}(r, t)   ←  comp.gas.Sigma             [g/cm²]  densidad real de gas
```

La fracción de cada especie acretada por el embrión se calcula directamente:

$$f_X^{\rm acc}(r,t) = \frac{\Sigma_{\rm dust,X}(r,t)}{\sum_Y \Sigma_{\rm dust,Y}(r,t)}$$

Los silicatos **no tienen su propio field** — se calculan como el remanente:

$$\Sigma_{\rm silicates}(r,t) = \max\!\left(0,\ \Sigma_{\rm dust,total}(r,t) - \sum_X \Sigma_{\rm dust,X}(r,t)\right)$$

**Método anterior (versión `SigmaIce_X`, ya retirado como primario):**

$$\Sigma_{\rm ice,X}(r,t) = f_X \cdot \Sigma_{\rm dust}(r,t) \cdot \theta\!\left[T(r,t) < T_{\rm sub,X}\right]$$

Este método era una estimación: usaba abundancias iniciales fijas y una función escalón en la temperatura, sin reflejar el drift y pile-up real de cada especie. El campo `SigmaDust_X` es el valor real calculado por tripodpy a través de la evolución de conservación de masa del componente.

**Especies soportadas:** H₂O, CO₂, CO (default). Extensible vía `pipeline.active_species`.

---

## 4. Por qué no se usaron las citas más recientes del v2_docs

El `PebbleAccretion_v2_docs.md` (generado por ChatGPT) incluía referencias más recientes: **McCloat+2025**, **Drążkowska+2022**, **Mulders+2021**, **Yap & Stevenson+2023**. No se incorporaron directamente en `PebbleAccretion2.py` por las siguientes razones:

### McCloat+2024/2025 ("PPOLs Model")
McCloat usa exactamente las mismas fórmulas de Ormel & Klahr (2010) y Lambrechts & Johansen (2012) como base. Su contribución es el análisis de arquitecturas planetarias y el mecanismo "pebble snow" (migración de snowlines que entrega hielo a planetas interiores), **no fórmulas nuevas para Ṁ_core o R_acc**. Citar McCloat en el código sin implementar el pebble-snow dinámicamente sería una cita falsa.

> **Pendiente de implementar:** el efecto pebble-snow (ya lo tenemos parcialmente, porque nuestro SigmaIce cambia con el tiempo, pero habría que estudiar el enriquecimiento de hielo en la zona interna cuando el snowline migra hacia adentro cruzando al embrión).

### Drążkowska+2022
Es un artículo de revisión (sin fórmulas nuevas propias). Las fórmulas que revisán son las mismas que ya citamos. Se puede añadir como contextualización general si el trabajo final se publica.

### Mulders+2021
Estudia el efecto de la masa estelar en el filtrado de pebbles por planetas gigantes. El código actual simula embriones aislados sin filtrado entre ellos. Incorporarlo requeriría simular múltiples embriones con reducción mutua del flujo.

### Yap & Stevenson+2023
Trata snowlines en discos **circum-planetarios** (alrededor de Júpiter/Saturno), no circum-estelares. Relevante para satélites, fuera del alcance actual.

### Efectos del v2_docs NO implementados (trabajo futuro)

| Efecto | Razón de exclusión | Prioridad |
|--------|-------------------|-----------|
| Excentricidad/inclinación en Δv | Sistema de embrión único, e≈0 válido | Media |
| Ablación en atmósferas | Requiere perfil atmosférico completo (Brouwers+2018) | Alta |
| Filtrado entre embriones | Requiere N embriones simultáneos | Alta |
| "Pebble snow" dinámico | Parcialmente implementado, mejorable | Alta |

---

## 5. Bugs corregidos respecto a versiones anteriores

| Bug | Versión afectada | Corrección |
|-----|-----------------|------------|
| `delta_v = eta * nu` (unidades incorrectas) | v1, ChatGPT | `delta_v = eta * v_K` |
| `grid/r` con shape `(Nt, Nr)` → error | v2 | `r_raw[0] if ndim==2 else r_raw` |
| Regímenes Bondi/Hill no separados | v1 | `min(R_H, b_Bondi, b_Hill)` |
| Transición 2D/3D brusca (if) | v1, ChatGPT | `erf(z)` continuo |
| `M_iso` sin correcciones α ni ∂lnP/∂lnr | v1, ChatGPT | Bitsch+2018 Eq. 5 completa |
| `enable_component_saving()` no guardaba Sigma | pipeline | `add_ice_sigma_fields()` con addfield() |
| Composición normalizada por bin grande (≈0 en r ext.) | v2 | Normalizar por Σ(SigmaIce_X) |
| Efecto atmósfera sin referencia y mal aplicado | ChatGPT | Eliminado, pendiente Brouwers+2018 |

---

## 6. Cómo reproducir los resultados

```bash
# 1. Correr el pipeline con SigmaIce guardado
python pipeline_snowlines.py
# → genera data_post_pipeline/pipeline_icefrac/ con 30 HDF5

# 2. Correr el módulo de acreción
python test_pebble_accretion.py
# → genera figs_pebbles/test_*.pdf

# 3. Ver tabla de composición
python diag_components.py
# → imprime tabla + guarda pebble_results.txt
```

---

## 7. Próximos pasos (Fase de Tesis)

1. **Simular hasta ~1 Myr** → obtener Waterworlds reales (f_H₂O > 10%)
2. **Variar masa estelar** (M★ = 0.5, 1.0, 1.5 M☉) → Fase 4 del plan de tesis
3. **Implementar filtrado** → cuando un embrión alcanza M_iso, reduce el flujo de pebbles para embriones más internos
4. **Diagrama composición vs. radio** → figura de publicación candidata
5. **Efecto pebble-snow** → cuantificar entrega de H₂O a la zona habitable cuando el snowline H₂O cruza al embrión

---

## 8. Referencias completas

| Referencia | Uso |
|------------|-----|
| Ormel & Klahr (2010), A&A 520, A43 | Radio de captura Bondi |
| Lambrechts & Johansen (2012), A&A 544, A32 | delta_v, radio Hill, flujo de pebbles |
| Youdin & Lithwick (2007), Icarus 192, 588 | H_peb settling |
| Liu & Ormel (2018), A&A 615, A178 | Transición 2D/3D erf |
| Lambrechts et al. (2014), A&A 572, A35 | Masa de aislamiento base |
| Bitsch et al. (2018), A&A 612, A30, Eq.5 | Correcciones f_α y f_P a M_iso |
| Haynes (1992); Fraser et al. (2001) | Abundancias moleculares iniciales |
| McCloat (2024, tesis) | Arquitecturas tipo-PPOLs, pebble-snow (contexto) |
| McCloat et al. (2025), arXiv:2509.14101 | Pebble snow, síntesis de población (contexto) |
| Drążkowska et al. (2022) | Revisión teoría de formación planetaria ALMA/Kepler (contexto) |
| Mulders et al. (2021) | Filtrado pebbles y masa estelar (futuro) |

---

## 9. Proceso actual con campos reales de componentes

Esta sección describe el flujo completo desde la simulación de disco hasta la acreción, usando los **campos reales** `SigmaDust_{sp}` y `SigmaGas_{sp}` introducidos en `add_ice_sigma_fields()`.

### 9.1 Paso 1 — Simulación del disco (`pipeline_snowlines.py`)

```python
pipeline = WaterworldPipeline("data_post_pipeline/mi_run")
pipeline.active_species = ["H2O", "CO2", "CO"]   # modificar según estudio
pipeline.setup_grid(rmin=1*c.au, rmax=300*c.au, Nr=200)
pipeline.setup_star()
pipeline.initialize_simulation()
pipeline.add_volatile_components()   # inyecta componentes desde chem.txt
pipeline.setup_physics()             # v_frag(T) dinámico por especie activa
pipeline.setup_star_evolution()      # R★ contrae de 2→1.5 R☉ en 10 Myr
pipeline.add_snowline_fields()       # grid/rsnow_{sp} en HDF5
pipeline.add_ice_sigma_fields()      # grid/SigmaDust_{sp} y SigmaGas_{sp} en HDF5
pipeline.sim.update()
pipeline.run_integration(t_end_years=1e5, num_snapshots=50)
```

**Resultado:** `Nt` archivos HDF5, cada uno con:
- `gas/Sigma`, `gas/T`, `gas/cs`, `gas/eta`, `gas/nu`
- `dust/Sigma` (Nr, Nbins), `dust/v/rad`, `dust/St`
- `grid/rsnow_H2O`, `grid/rsnow_CO2`, `grid/rsnow_CO`
- `grid/SigmaDust_H2O`, `grid/SigmaDust_CO2`, `grid/SigmaDust_CO`
- `grid/SigmaGas_H2O`, `grid/SigmaGas_CO2`, `grid/SigmaGas_CO`

### 9.2 Paso 2 — Acreción de pebbles (`PebbleAccretion2.py`)

```python
from PebbleAccretion2 import PebbleAccretionModule

pam = PebbleAccretionModule.from_datadir("data_post_pipeline/mi_run", M_star=1.0)
# → detecta: Especies detectadas: ['H2O', 'CO2', 'CO'] + silicates
# → silicates = dust_total - SigmaDust_H2O - SigmaDust_CO2 - SigmaDust_CO

results = pam.run_growth([1.0, 2.0, 3.0, 5.0, 8.0, 15.0], M0_g=7.3e25)
pam.summary(results)
```

### 9.3 Diferencia clave respecto a la versión anterior

| Aspecto | Método antiguo (`SigmaIce_X`) | Método actual (`SigmaDust_X`) |
|---|---|---|
| Origen | Estimado: `f_X × Σ_dust × θ(T<T_sub)` | Real: `comp.dust.Sigma.sum(-1)` de tripodpy |
| Drift por especie | Ignorado | Incluido (tripodpy integra masa por componente) |
| Pile-up en snowlines | Solo vía θ escalón en T | Dinámico: emerge de condensación/evaporación |
| Silicatos | Fracción fija `f_sil = 2e-3/total` | Remanente: `Σ_dust_total − Σ_volátiles` |
| Especies detectadas | Hardcodeadas: H₂O, CO₂, silicatos | Auto-detectadas desde claves HDF5 |

### 9.4 Perfil `v_frag(T)` por especie activa

`setup_physics()` construye el perfil automáticamente a partir de `active_species`:

| Temperatura | Especie | v_frag | Referencia |
|---|---|---|---|
| T ≥ 150 K | Silicatos (baseline, siempre) | 100 cm/s | Birnstiel+2012 |
| T < 150 K | H₂O ice (si activa) | 1000 cm/s | Gundlach & Blum 2015 |
| T < 70 K | CO₂ ice (si activa) | 500 cm/s | Gundlach+2018 |
| T < 25 K | CO ice (si activa) | 300 cm/s | Dominik & Tielens 1997 |

Si se elimina una especie de `active_species`, esa banda de temperatura cae al baseline de silicatos (100 cm/s).

El algoritmo aplica los escalones de mayor a menor T_sub (H₂O → CO₂ → CO), de modo que cada `np.where` subsiguiente solo sobreescribe la zona más fría, sin afectar la zona exterior.

### 9.5 Parámetros configurables en `WaterworldPipeline`

| Atributo | Descripción | Default |
|---|---|---|
| `active_species` | Lista de especies activas (deben estar en chem.txt) | `["H2O","CO2","CO"]` |
| `vfrag_params` | `{sp: (T_sub_K, v_frag_cm_s)}` para cada especie | H₂O/CO₂/CO estándar |
| `vfrag_silicates` | Baseline de silicatos [cm/s] | `100.0` |
| `M_star_ini` | Masa estelar inicial | `1.0 M☉` |
| `R_star_ini` | Radio estelar inicial (joven, inflado) | `2.0 R☉` |
| `T_star_ini` | Temperatura efectiva | `5778 K` |

### 9.6 Fix del integrador de Euler (`M_iso` clamp)

El integrador de `run_growth` usa Euler explícito con los snapshots como pasos de tiempo. Sin control, un `dt` grande podía hacer que `M_core` excediera `M_iso` en un solo paso (overshoot). El fix aplica:

```python
dM = min(Mdot_core * dt, max(0.0, M_iso - M_core))
```

Esto garantiza que `M_core` nunca supere `M_iso` independientemente de cuántos snapshots se usen.

> **Nota:** más snapshots → `dt` más pequeño → integración más precisa (menos overshoot en masa). El `M_iso` clamp es una salvaguarda adicional.

---

## 10. Visualización: `plot_diagnostics.py` (actualizado)

### 10.1 Qué cambió

| Método | Antes | Ahora |
|---|---|---|
| `plot_gas_components` | Requería `frame.dmp` (solo estado final) | Lee `grid/SigmaGas_X` del HDF5 elegido con `it` |
| `plot_dust_components` | Requería `frame.dmp` | Lee `grid/SigmaDust_X` del HDF5; silicatos = remanente |
| `_detect_components` | Solo escaneaba `components.*` del dump | También detecta `SigmaDust_X` en `grid/` → `self.volatile_species` |
| `plot_hovmoller_comp` | No existía | Nuevo: Hovmöller por especie sobre tiempo |
| `__main__` | 7 outputs | 8 outputs (añade `hovmoller_comp_dust` y `hovmoller_comp_gas`) |

### 10.2 Uso

```python
from plot_diagnostics import SnowlineDiagnostics

d = SnowlineDiagnostics(
    "data_post_pipeline/pipeline_v3_Sigma_update",
    savedir="figs",
    r_trim=0.93,
    t_min_yr=100.0
)

# Al inicializar imprime automáticamente:
# → Especies con Sigma en HDF5: ['H2O', 'CO2', 'CO']

# Snapshot específico (ya no necesita frame.dmp):
d.plot_gas_components(it=-1)    # último snapshot
d.plot_gas_components(it=0)     # snapshot inicial
d.plot_dust_components(it=15)   # snapshot intermedio

# Hovmöller por especie — evolución temporal completa:
d.plot_hovmoller_comp(quantity='dust')   # SigmaDust_H2O / CO2 / CO / silicatos
d.plot_hovmoller_comp(quantity='gas')    # SigmaGas_H2O / CO2 / CO
```

### 10.3 Outputs generados por `__main__`

| N° | Archivo PDF | Descripción |
|---|---|---|
| 1 | `hovmoller_eps.pdf` | ε = Σ_dust/Σ_gas (espacio-tiempo) |
| 2 | `hovmoller_Sigma_dust.pdf` | Σ_dust total (espacio-tiempo) |
| 3 | `size_distribution_it-001.pdf` | Distribución de tamaño (último snapshot) |
| 4 | `pebble_flux_hovmoller.pdf` | Ṁ_pebble(r, t) en M⊕/yr |
| 5 | `profiles_it-001.pdf` | η, St, a_max, Σ gas+dust (4 paneles) |
| 6a | `gas_components_it-001.pdf` | Σ_gas por especie (H₂O, CO₂, CO + total) |
| 6b | `dust_components_it-001.pdf` | Σ_dust por especie + silicatos remanentes |
| 7a | `hovmoller_comp_dust.pdf` | Hovmöller de Σ_dust: panel por especie + silicatos |
| 7b | `hovmoller_comp_gas.pdf` | Hovmöller de Σ_gas: panel por especie |
| 8 | `hovmoller_rt.pdf` | Σ_gas, T, a_max multi-panel R(t) |

### 10.4 `plot_hovmoller_comp` — descripción del output

Para `quantity='dust'` con `active_species = ['H2O', 'CO2', 'CO']` genera **4 paneles apilados**:

```
Panel 1: SigmaDust_H2O(r, t)   — cmap magma
Panel 2: SigmaDust_CO2(r, t)   — cmap magma
Panel 3: SigmaDust_CO(r, t)    — cmap magma
Panel 4: Silicatos(r, t)       — cmap cividis
         = dust_total - SigmaDust_H2O - SigmaDust_CO2 - SigmaDust_CO
```

- Eje X: tiempo [años, escala log]
- Eje Y: radio [AU, escala log]
- Snowlines superpuestas en cada panel (líneas `--` coloreadas)
- Título del panel en el color de la especie (COMP_COLORS)

### 10.5 Mapeo correcto de índices temporales

El método `_resolve_it(it, Nt)` devuelve un índice en el subconjunto filtrado por `t_min_yr`. Para acceder al HDF5 completo (shape `(Nt_full, Nr)`), se usa:

```python
it_full = int(np.where(self.t_mask)[0][it_abs])
sig = np.asarray(field_data)[it_full, self.r_mask]
```

Esto evita el bug anterior donde `it_abs` se usaba directamente en el array completo cuando `t_min_yr > 0` filtraba snapshots del inicio.

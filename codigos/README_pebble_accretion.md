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
  └── add_ice_sigma_fields()          ← guarda grid/SigmaIce_X en cada HDF5
  └── add_snowline_fields()           ← guarda grid/rsnow_X en cada HDF5
  └── run_integration()               ← genera data_post_pipeline/pipeline_icefrac/

PebbleAccretion2.py
  └── PebbleAccretionModule
        ├── from_datadir(datadir)     ← apila 30 HDF5 en arrays (Nt, Nr, ...)
        │     Detecta grid/SigmaIce_X (nuevo) o components/*/dust/Sigma (legacy)
        ├── run_growth(embryos_AU)    ← integra dM/dt para cada embrión
        │     Composición: normaliza por sum(SigmaIce_X) [partición de unidad]
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

### 3.6 Composición química

Campo `grid/SigmaIce_X` guardado en cada snapshot HDF5:

$$\Sigma_{\rm ice,X}(r,t) = f_X \cdot \Sigma_{\rm dust}(r,t) \cdot \theta\!\left[T(r,t) < T_{\rm sub,X}\right]$$

Fracción acretada por cada especie:

$$f_X^{\rm acc} = \frac{\Sigma_{\rm ice,X}}{\sum_Y \Sigma_{\rm ice,Y}}$$

Abundancias iniciales (Haynes 1992; Fraser et al. 2001): H₂O: 1.6×10⁻⁴, CO₂: 4×10⁻⁵, CO: 8×10⁻⁵, silicatos: 2×10⁻³

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

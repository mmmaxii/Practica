# `PebbleAccretion2.py` — Documentación Física  `#CORREGIDO`

> Este documento fue revisado y corregido respecto a la versión original generada por ChatGPT.  
> Las correcciones se marcan con **`[CORREGIDO]`** y están justificadas con la referencia bibliográfica exacta.

---

## 1. Velocidad relativa pebble-embrión Δv

**`[CORREGIDO]`** La versión anterior usaba `delta_v = eta * nu`, que es **dimensionalmente incorrecto** (`nu` tiene unidades cm²/s, no cm/s).

La fórmula correcta es:

$$\Delta v = \eta \, v_K, \qquad v_K = \Omega_K \cdot r$$

donde $\eta$ es el parámetro adimensional de gradiente de presión del gas (disponible directamente como `gas.eta` en tripodpy).

**Referencia**: Lambrechts & Johansen (2012), A&A 544, A32, **Eq. 16**.

---

## 2. Radio de captura efectivo $R_{\rm acc}$

**`[CORREGIDO]`** La versión anterior no diferenciaba entre los regímenes de Bondi y Hill. La implementación correcta usa ambos y toma el mínimo:

$$R_{\rm acc} = \min\left(R_H,\; b_{\rm Bondi},\; b_{\rm Hill}\right)$$

donde:

$$R_H = r\left(\frac{M}{3M_*}\right)^{1/3} \quad \text{(radio de Hill, estándar)}$$

$$b_{\rm Bondi} = \sqrt{St} \cdot \frac{G\,M}{\Delta v^2}$$
> **Referencia**: Ormel & Klahr (2010), A&A 520, A43, **Eq. A9** (régimen de Bondi)

$$b_{\rm Hill} = \left(\frac{St}{0.1}\right)^{1/3} R_H$$
> **Referencia**: Lambrechts & Johansen (2012), A&A 544, A32, **Eq. A10** (régimen de Hill)

El mínimo garantiza que no se use una escala mayor que la esfera de influencia gravitacional real.

---

## 3. Altura del disco de pebbles $H_{\rm peb}$

**`[CORREGIDO]`** La versión anterior usaba `self.H` (array global) como escalar en una comparación condicional → comportamiento indefinido.

La fórmula correcta, interpolada en $r_{\rm emb}$:

$$H_{\rm peb} = H_{\rm gas}\,\sqrt{\frac{\alpha}{\alpha + St}}$$

**Referencia**: Youdin & Lithwick (2007), Icarus 192, 588, **Eq. 12**.

Donde $\alpha = \nu/(c_s H_{\rm gas})$, estimado directamente de los campos `gas.nu`, `gas.cs` del HDF5 de tripodpy.

---

## 4. Tasas de acreción 2D y 3D — transición suave

La transición brusca `if R_acc > H_peb` fue **`[CORREGIDO]`** por una interpolación continua con función error:

$$\dot{M}_{\rm 2D} = 2\,R_{\rm acc}\,\Sigma_{\rm peb}\,\Delta v$$

$$\dot{M}_{\rm 3D} = \pi\,R_{\rm acc}^2\,\rho_{\rm peb}\,\Delta v, \qquad \rho_{\rm peb} = \frac{\Sigma_{\rm peb}}{\sqrt{2\pi}\,H_{\rm peb}}$$

$$\dot{M}_{\rm core} = \mathrm{erf}(z)\,\dot{M}_{\rm 2D} + \left[1 - \mathrm{erf}(z)\right]\,\dot{M}_{\rm 3D}, \qquad z = \frac{R_{\rm acc}}{\sqrt{2}\,H_{\rm peb}}$$

**Referencia**: Liu & Ormel (2018), A&A 615, A178, **Appendix A.1**.

> Nota: la versión anterior también omitía el factor de eficiencia de captura. La fórmula geométrica $2 R_{\rm acc} \Sigma_{\rm peb} \Delta v$ con el $R_{\rm acc}$ correcto del settling regime es la expresión estándar de la literatura y no requiere factor adicional cuando $R_{\rm acc}$ se calcula correctamente (O&K 2010; L&J 2012).

---

## 5. Masa de aislamiento $M_{\rm iso}$

**`[CORREGIDO]`** La versión anterior usaba solo la fórmula simplificada de Lambrechts+2014. La implementación correcta usa Bitsch et al. (2018) con correcciones de viscosidad y gradiente de presión:

$$M_{\rm iso} = 25\,M_\oplus\left(\frac{h}{0.05}\right)^3\,f_\alpha\,f_P$$

donde:

$$f_\alpha = 0.34\left[\log_{10}\!\left(\frac{\alpha}{10^{-3}}\right) + 3\right] + 0.66$$

$$f_P = 1 - \frac{\partial\ln P/\partial\ln r + 2.5}{6}, \qquad \frac{\partial\ln P}{\partial\ln r} = \frac{-2\eta}{(H/r)^2}$$

La relación $\partial\ln P/\partial\ln r = -2\eta/(H/r)^2$ se deriva de la definición de $\eta$ y permite usar directamente `gas.eta` y `gas.cs` del HDF5.

**Referencias**:
- Fórmula base: Lambrechts et al. (2014), A&A 572, A35
- Correcciones completas: **Bitsch et al. (2018), A&A 612, A30, Eq. 5**

Como cota inferior conservadora se mantiene $M_{\rm iso} \geq 20\,M_\oplus\,(h/0.05)^3$.

---

## 6. Lectura del grid HDF5

**`[CORREGIDO]`** `grid/r` en el HDF5 de tripodpy puede tener shape `(Nt, Nr)` al ser leído con `read.all()`. La versión anterior usaba directamente `f['grid/r'][:]` sin manejar ambas posibilidades, lo que producía que `Omega_K` tuviera shape `(Nt, Nr)` y rompía toda la interpolación subsecuente.

Corrección:
```python
r_raw  = f['grid/r'][:]
self.r = r_raw[0] if r_raw.ndim == 2 else r_raw   # siempre (Nr,)
```

---

## 7. Efecto atmósfera — eliminado

**`[CORREGIDO]`** La versión anterior aplicaba `Mdot *= f_silicates` cuando `M_core > 0.1 M_⊕`, bajo el argumento de que el hielo se evapora en la atmósfera. Esto fue **eliminado** por las siguientes razones:
1. El umbral 0.1 M_⊕ es arbitrario y no tiene referencia bibliográfica sólida.
2. El efecto se aplicaba al `Mdot` total en vez de solo a la fracción de hielo, reduciendo también la acreción de silicatos.
3. El seguimiento de composición ya refleja naturalmente los cambios cuando el embrión cruza un snowline (la fracción $f_X$ de cada especie cambia según `components.X.dust.Sigma`).

Si se desea modelar sublimación en atmósfera, debe implementarse como un post-proceso separado con física de la ablación (ej. Brouwers & Ormel 2020).

---

## 8. Flujo de pebbles: velocidad ponderada

El flujo $\dot{M}_{\rm peb} = 2\pi r \Sigma_{\rm peb} |v_r|$ usa la velocidad del mayor Stokes (`vr[:, -1]`), correspondiente al grano de mayor drift. Esto es una simplificación conservadora: en un modelo de dos poblaciones, el bin grande domina el transporte de masa.

Una mejora futura sería usar la velocidad promedio ponderada por masa:
$$v_r^{\rm eff} = \frac{\sum_k \Sigma_k v_{r,k}}{\sum_k \Sigma_k}$$

**Referencia del flujo**: Lambrechts & Johansen (2012), A&A 544, A32.

---

## 9. Cruce de snowlines — implementación automática

El modelo no requiere manejo explícito del cruce: las fracciones $f_X(r_{\rm emb}, t) = \Sigma_X / \Sigma_{\rm peb}$ calculadas desde los campos `components.X.dust.Sigma` del HDF5 ya son 0 cuando el embrión está dentro del snowline (el hielo se ha sublimado y el gas del componente no contribuye a `dust.Sigma`). Cuando la snowline migra y el embrión queda fuera, $f_X > 0$ automáticamente.

Para diagnóstico, se registra `rsnow_H2O(t)` en el histórico de cada embrión (columna 6 del array de salida).

---

## Referencias completas

| Referencia | Uso en el código |
|---|---|
| Ormel & Klahr (2010), A&A 520, A43 | Radio de captura Bondi: $b_{\rm Bondi} = \sqrt{St}\cdot GM/\Delta v^2$ |
| Lambrechts & Johansen (2012), A&A 544, A32 | $\Delta v = \eta v_K$; radio Hill: $b_{\rm Hill} = (St/0.1)^{1/3}R_H$; flujo de pebbles |
| Youdin & Lithwick (2007), Icarus 192, 588 | Altura del disco de pebbles: $H_{\rm peb} = H_{\rm gas}\sqrt{\alpha/(\alpha+St)}$ |
| Liu & Ormel (2018), A&A 615, A178 | Transición suave 2D/3D via $\rm erf$ |
| Lambrechts et al. (2014), A&A 572, A35 | Masa de aislamiento base |
| Bitsch et al. (2018), A&A 612, A30, Eq. 5 | Correcciones $f_\alpha$, $f_P$ a $M_{\rm iso}$ |
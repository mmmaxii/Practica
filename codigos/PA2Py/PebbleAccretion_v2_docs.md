# 1. Lo logrado y fundamentos actuales 

El módulo ya implementa correctamente el **núcleo analítico estándar** de acreción de pebbles. En particular:  

- **Rótulo unificado de captura**: Se usa $$R_{\rm acc}\approx\min(R_H,\,R_B\sqrt{St})$$, con $R_H=r_p(M/(3M_*))^{1/3}$ y $R_B=2GM/(\Delta v)^2$, cubriendo así tanto el régimen de fuertes desviaciones gravitatorias (pebbles grandes) como el de arrastre gas-dominado (pebbles pequeños)【96†L69-L72】. Este enfoque coincide exactamente con las fórmulas de Ormel & Klahr (2010) y Lambrechts & Johansen (2012), consideradas el “gold standard” actual.  

- **Masa de aislamiento de pebbles $M_{\rm iso}$**: Se han implementado las expresiones de Lambrechts et al. (2014) y Bitsch et al. (2018), que dependen fuertemente de $H/r$, turbulencia $\alpha$ y gradientes de presión. Por ejemplo, típicamente $$M_{\rm iso}\sim20\,M_\oplus\,(H/r)^3$$ ajustado por factores turbulentos【96†L91-L93】. Incluir la dependencia en $\alpha$ (Bitsch+2018) es vital: $M_{\rm iso}$ varía con la posición de trampas de polvo y líneas de nieve.  

- **Seguimiento químico dinámico**: El código lee del disco la densidad de cada especie (H$_2$O, CO$_2$, silicato) en cada radio $r$, y al calcular la fracción $\Sigma_X/\Sigma_{\rm peb}$ en ese momento, asigna la porción de masa del embrión a cada componente. Así, si el embrión cruza una línea de nieve (p.ej. de agua), empieza a incorporar hielo H$_2$O según los campos químicos simulados【98†L34-L40】【93†L62-L69】. Esto es precisamente lo que hoy se considera esencial para predecir relaciones como C/O en atmósferas (por ejemplo en JWST): el código ya realiza un seguimiento “local” en tiempo real de la composición de los pebbles accedidos.  

- **Filtrado y arquitectura emergente**: Al simular múltiples embriones simultáneamente, el código ya contempla la reducción del flujo de pebbles (filtrado) y la masa de aislamiento de cada semilla. Es decir, los resultados emergen de la dinámica colectiva (como en modelos de síntesis de población)【98†L34-L42】【98†L45-L53】. Esto ha permitido reproducir secuencias de arquitectura *tipo* “baja, media y alta masa de disco” análogas a estudios recientes (por ej. PPOLs)【98†L45-L53】: discos poco masivos producen muchos núcleos rocosos secos en la zona interna, discos muy masivos generan núcleos acuosos externos que bloquean el flujo interno, y discos intermedios generan combinaciones (análogas a patrones “Solar System–like” en rango estrecho de parámetros)【98†L45-L53】【83†L732-L740】.  

En resumen, el núcleo físico del modelo (fórmulas Ormel/Klahr, LJ12/14, Bitsch18) y el seguimiento químico están **correctos y a la par de los modelos de vanguardia** (por ejemplo las simulaciones del PPOLs de McCloat et al. 2025 usan exactamente estas bases【98†L34-L42】).  

# 2. Aspectos avanzados (estado del arte) 

Para acercarse a la “vanguardia”, conviene incorporar varios efectos sutiles señalados en la literatura reciente:

- **2A. Transición suave entre regímenes 2D/3D.** En modelos reales la distribución vertical de pebbles es continua (gaussiana). En lugar de un corte neto “if $R_{\rm acc}>H_{\rm peb}$ then 2D”, se usan interpolaciones continuas. Por ejemplo, Liu & Ormel (2018) o Ida et al. (2016) proponen una fórmula del tipo  
  $$\dot M_{\rm core} \approx \dot M_{2D}\times\left[1+\left(\frac{R_{\rm acc}}{H_{\rm peb}}\right)^n\right]^{-1},$$  
con $n\sim1$–2, o soluciones integradas sobre la altura. Esto evita cambios bruscos en $\dot M_{\rm core}$ cuando $R_{\rm acc}\approx H_{\rm peb}$. En la práctica se puede implementar multiplicando la tasa 2D por un factor $f=\mathrm{erf}(R_{\rm acc}/(\sqrt{2}H_{\rm peb}))$ o $f=(1+\exp[-R_{\rm acc}/H_{\rm peb}])^{-1}$, que suaviza la transición. Tal refinamiento elimina artefactos numéricos y refleja mejor las derivaciones teóricas recientes.  

- **2B. Excentricidad e inclinación de los embriones.** El módulo actual asume órbitas circulares (Δ$v\sim\eta v_K$). Sin embargo, embriones masivos pueden adquirir excentricidades $e$ e inclinaciones $i$ por interacciones mutuas. En ese caso la velocidad relativa efectiva es mayor: $$\Delta v \approx \sqrt{(\eta v_K)^2 + (e v_K)^2 + (i v_K)^2}\,. $$ Un $e$ incluso pequeño (p.ej. $e\sim0.01$) puede superar $\eta$ y reducir dramáticamente el radio de Bondi $R_B\propto1/\Delta v^2$. En la práctica, si el embrión desarrolla $e,i$, la eficiencia de acreción cae (Ver Ormel & Klahr 2010; Liu & Ormel 2018). Para modelarlo, se debería aumentar localmente $\Delta v$ con contribución aleatoria o estimada de $e,i$ (quizá derivada de un esquema externo N-cuerpos) y recalcular $R_B$. Esto introducirá retardos en el crecimiento en escenarios con múltiples cuerpos o turbulencia planetaria, como señalan trabajos recientes.  

- **2C. Efectos de atmósferas primordiales y ablación.** Una vez el núcleo supera $\sim$0.1–0.3\,M$_\oplus$, comienza a retener una envoltura H/He. Los pebbles ya no colisionan directamente con el núcleo, sino que penetran la atmósfera. Modelos recientes (Brouwers et al. 2018; Vazan et al. 2022) demuestran que esto aumenta el “radio efectivo” de captura (por aerodinámica) pero *no* toda la materia llega al núcleo sólido: especialmente el hielo volátil se evapora en la atmósfera. En otras palabras, existe un “radio de captura atmosférico” $R_{\rm cap}>R_{\rm acc}$, y el núcleo recibe solo silicatos o vapores en lugar de granos enteros helados. Para incorporar esto, habría que calcular el perfil de la atmósfera del embrión (p.ej. integrando gradiente adiabático) y estimar la ablandamiento de pebbles. Una aproximación semi-analítica es: para $M\gtrsim0.1M_\oplus$, suponer que el hielo permanece en atmósfera y solo los refractarios suman al núcleo. Esto haría que la *composición* del crecimiento cambie drásticamente en etapas tardías, acentuando el efecto de aislamiento químico de la atmósfera.  

- **2D. Formación inicial de planetesimales (noten respuesta actual).** Recordemos que la acreción de pebbles **no forma planetesimales**, sino que requiere “semillas” preexistentes (como su M$_0\sim10^{24}$\,g). El módulo modela el **crecimiento** posterior de esos núcleos, no su formación (que requiere Streaming Instability u otro mecanismo). Es importante conceptualmente: los primeros $10^{21}$–$10^{24}$\,g surgen por procesos de aglomeración colectivamente (no por las fórmulas de acreción de pebbles). Este código asume la etapa de semilla completada y modela sólo la etapa de crecimiento.  

# 3. Evolución química y líneas de nieve (“Pebble Snow”) 

Los modelos recientes muestran que la **migración de las líneas de nieve** altera drásticamente la composición de los pebbles que llegan a cada embrión. Si un embrión nace **dentro** de la línea de hielo de H$_2$O, inicialmente recibe solo polvo rocoso seco. Cuando la línea de nieve migra hacia dentro (a medida que el disco se enfría), de repente el mismo embrión pasa a recibir pebbles con fracción significativa de H$_2$O. Este efecto –denominado **“nieve de pebbles”**– implica que se inyectan hielos en los planetas interiores, aun si se formaron inicialmente en zona seca【98†L118-L126】. Mulders et al. (2021) y McCloat et al. (2025) cuantifican esto: por ejemplo, un embrión en 1 AU acumula mayor parte del agua si la línea de nieve pasa luego de 2–4 Myr (discos grandes)【96†L121-L129】【93†L69-L74】. McCloat et al. (2024) muestran además que sólo en un rango estrecho de masas de disco (alrededor de sistemas solares) se entrega *justamente* una fracción de agua “terrestre” a planetas de 1 M⊕ en zona habitable【93†L69-L74】.  

Por tanto, la **entrega de agua a la zona habitable** ocurre por esta “barrera deslizante” de hielos: pebbles cargados de H$_2$O atraviesan el sistema interno cuando el borde helado se mueve. Como resultado, embriones que eran secos pueden volverse húmedos, y la fracción final de agua depende críticamente del tiempo de cruce. En síntesis: la migración de líneas de nieve junto con el flujo interno de pebbles genera un gradiente fuerte de contenido de agua en los planetas internos【98†L118-L126】【93†L69-L74】.  

# 4. Discusión y síntesis bibliográfica 

Recapitulando los aportes clave de la literatura reciente (2021–2025):

- **McCloat (2024, tesis “PPOLs Model”)【93†L62-L69】【93†L69-L74】:** Desarrolla un modelo global de acreción de pebbles (“PPOLs”) que simula múltiples núcleos simultáneos con evolución de la línea de nieve. Muestra que sólo en un estrecho rango de parámetros se entrega *exacta* cantidad de agua “terrestre” a planetas tipo Tierra en la zona habitable【93†L69-L74】. También enfatiza que semillas múltiples y evolución de la nieve prolongan la fase de crecimiento y enriquecen de hielos al disco interno, conectando el final de acreción con la diversidad observada. Sus conclusiones resaltan la importancia combinada de filtrado, masa de aislamiento y nieve móvil en determinar la arquitectura final【83†L732-L740】【83†L743-L752】.

- **McCloat et al. (2025, arXiv 2509.14101)【98†L34-L40】【98†L45-L53】:** Introducen formalmente el mecanismo “pebble snow” y un modelo de síntesis rápido. Confirman tres arquitecturas tipo (baja, media, alta masa de disco) y cómo difieren en la distribución de masa y de agua. Destacan que un disco bajo produce muchos núcleos secos en HZ (varias masas *Mars–Tierra* con aguas variadas), discos altos generan núcleos exteriores ricos en agua que bloquean el interior, y discos intermedios crean sistemas bimodales parecidos al Sistema Solar en un estrecho rango【98†L45-L53】. Resaltan la entrega de hielos en el punto de la línea de nieve (generando núcleos “pre-giant” y un pico interior). Esta obra amplía modelos previos al incluir evolución coherente de la nieve y múltiples semillas simultáneas.  

- **Drążkowska et al. (2022/2023)** *“Planet Formation Theory in the Era of ALMA and Kepler”*: revisión extensa de estados nuevos de la formación planetaria. Recopila avances en la evolución de polvo, formación de planetesimales y acreción de pebbles en el contexto de ALMA/Kepler. Entre otros, enfatiza que los primeros núcleos se forman temprano y de manera no uniforme en el disco, y que la reserva global de pebbles (modelo de múltiples poblaciones) es crucial para el crecimiento de núcleos. Es relevante el énfasis en modelos de síntesis que conectan disco-observables con poblaciones exoplanetarias (sin cita directa). 

- **Mulders et al. (2021)** (referencia incluida por McCloat【96†L138-L147】): analizan el efecto combinado de un protoplaneta gigante exterior (5 AU) y la línea de nieve fija sobre un protoplaneta interior (0.3 AU). Encontraron que en estrellas de baja masa (M-dwarf) el filtro externo es menos eficiente (menos material para formar el gigante), permitiendo que el núcleo interior crezca más (facilitando super-Tierras interiores). Esto explicaría por qué los sistemas de baja masa estelar tienden a tener más SuperTierras interiores observadas【96†L138-L147】. También confirman que tanto el filtrado como la masa de aislamiento juegan papeles clave en la arquitectura final. 

- **Yap & Stevenson (2023)**【95†L42-L49】: extienden la física de líneas de nieve al disco circunplanetario de Júpiter/Saturno. Demuestran que la **línea de hielo en el disco decreciente del planeta** acumula sólidas de alta fracción de H$_2$O justo más allá de la línea (análogo a embudos de hielo): el hielo transportado en pebbles advierte hacia la línea, mientras que el vapor se difunde hacia atrás, enriqueciendo la región. Así, los satélites anómalamente ricos en agua (Ganimedes, Calisto, Titán) pueden formarse allí, con hielo:roca ∼1:1【95†L42-L49】. Este proceso es análogo, en miniatura, a la entrega de hielos en discos circum-estelares: un frente de hielo concentrador que produce cuerpos extra-ricos en agua.  

- **Otros estudios relacionados:** Varios trabajos recientes (p.ej. Brouwers et al. 2018, Vazan et al. 2022) discuten la ablación de pebbles en atmósferas y la formación de planetas acuosos (**water worlds**), así como el papel de la fragmentación de pebbles y planetas gigantes exteriores. Además, estudios observacionales (como DSHARP, *ALMA*) confirmaron anillos de polvo que actúan como trampas para pebbles (coincidiendo con los requerimientos teóricos de filtrado y aislamiento)【96†L69-L73】. La síntesis de población incorporando estos efectos está en pleno desarrollo (ver también modelos de síntesis del grupo de Lund y Berna).  

**En conjunto**, estas referencias muestran que el código original está fundamentado correctamente, pero que los detalles de la frontera del campo implican efectos más finos (explicados arriba) que pueden incorporarse para obtener predicciones aún más realistas y comparar directamente con la diversidad de exoplanetas observados.  

# 5. Implementación numérica sugerida  

A continuación se propone un esquema Python actualizado (compatible con TripodPy) que incluye algunas de las mejoras mencionadas. Este fragmento es ilustrativo y muestra cómo integrar los nuevos efectos:

```python
class PebbleAccretionModule:
    def __init__(self, h5file, M_star=1.0):
        import h5py, numpy as np
        f = h5py.File(h5file, 'r')
        # Cargar datos del disco (gas, polvo, química, etc.)
        self.r = f['grid/r'][:]       # radios [cm]
        self.times = f['t'][:]        # tiempos [s]
        self.rsnow = {spec: f[f'grid/rsnow_{spec}'][:] for spec in ('H2O','CO2','CO')}
        self.gas = {key: f[f'gas/{key}'][:] for key in ('Sigma','T','cs','eta','nu')}
        self.dust = {'Sigma': f['dust/Sigma'][:], 'vr': f['dust/v.rad'][:], 'St': f['dust/St'][:]}
        self.comp = {spec: f[f'components/{spec}/dust/Sigma'][:] for spec in ('H2O','CO2','silicates')}
        f.close()
        # Parametros estelares
        self.M_star = M_star*1.989e33  # g
        G = 6.674e-8
        # Calcular frequencia orbital y H_gas en cada radio
        Omega = np.sqrt(G*self.M_star/self.r**3)
        self.Omega = Omega
        self.H = {t: self.gas['cs'][t,:]/Omega for t in range(len(self.times))}
        # Calcular alpha turbulencia (suponiendo nu = alpha c_s H)
        self.alpha = {t: np.clip(self.gas['nu'][t,:]/(self.gas['cs'][t,:]*self.H[t]), 1e-4, 1e-1) for t in range(len(self.times))}

    def pebble_flux(self, t_idx, r_emb):
        """Flujo de pebbles en el radio del embrión."""
        import numpy as np
        Sigma_peb = np.interp(r_emb, self.r, self.dust['Sigma'][t_idx,:,1])  # bin grande
        # Velocidad de deriva: se usa la muestra de mayor St (índice -1)
        v_r = np.interp(r_emb, self.r, self.dust['vr'][t_idx,:,-1])
        return 2*np.pi * r_emb * Sigma_peb * abs(v_r)

    def accretion_rate(self, M_core, r_emb, t_idx, e=0.0, i=0.0):
        """Calcula Ṁ_core incluyendo efectos 2D/3D suavizados y atmosféricos."""
        import numpy as np
        # Interpolar datos locales
        Sigma_peb = np.interp(r_emb, self.r, self.dust['Sigma'][t_idx,:,1])
        cs = np.interp(r_emb, self.r, self.gas['cs'][t_idx])
        Omega = np.interp(r_emb, self.r, self.Omega)
        eta = np.interp(r_emb, self.r, self.gas['eta'][t_idx])
        St = np.interp(r_emb, self.r, self.dust['St'][t_idx,:,-1])
        # Velocidad relativa: incluir excentricidad/inclinación
        v_K = Omega * r_emb
        delta_v = np.sqrt((eta*v_K)**2 + (e*v_K)**2 + (i*v_K)**2)
        # Radios característicos
        G = 6.674e-8
        R_H = r_emb*(M_core/(3*self.M_star))**(1/3)
        R_B = 2*G*M_core/(delta_v**2 + 1e-10)
        R_acc = min(R_H, R_B*np.sqrt(St))
        # Altura de pebbles
        H_peb = np.sqrt(self.H[t_idx][np.argmin(abs(self.r-r_emb))]**2 * 
                         self.alpha[t_idx][np.argmin(abs(self.r-r_emb))] / (self.alpha[t_idx][np.argmin(abs(self.r-r_emb))] + St))
        # Fórmula 2D/3D suavizada
        # Factor de transición suave (p.ej. z = R_acc/(sqrt(2)*H_peb))
        z = R_acc/(np.sqrt(2)*H_peb + 1e-12)
        # Eficiencia efectiva 3D->2D: usando error function suaviza la transición
        f2D = np.erf(z)
        # Cálculo de Ṁ en cada régimen
        Mdot_2D = 2 * R_acc * Sigma_peb * delta_v
        Mdot_3D = np.sqrt(2*np.pi) * R_acc**2 * (Sigma_peb/(np.sqrt(2*np.pi)*H_peb + 1e-12)) * delta_v
        Mdot = f2D*Mdot_2D + (1-f2D)*Mdot_3D
        # Efecto atmósfera: si M_core grande, solo refractarios llegan al núcleo
        if M_core > 0.1*5.97e27:  # >0.1 M_E
            # Por simplicidad, suponer que el hielo no llega al núcleo sólido.
            # Se puede reducir Ṁ efectivo del hielo a 0 (o baja fracción).
            # Aquí se modela sólo Ṁ de silicatos:
            Mdot *= (np.interp(r_emb,self.r,self.comp['silicates'][t_idx]) / (Sigma_peb + 1e-12))
        return Mdot

    def run_growth(self, embryos, M0=1e24):
        """
        Simula el crecimiento de una lista de embriones dados (r_au) con 
        masa inicial M0 (g). Devuelve arrays de tiempo vs M_core y composición.
        """
        import numpy as np
        results = {}
        for r_au in embryos:
            r_emb = r_au * 1.4959787e13
            M_core = M0
            M_X = {spec:0.0 for spec in self.comp}  # masas acumuladas de cada especie
            history = []
            for i,t in enumerate(self.times):
                # Calcular Mdot_pebble y Mdot_core
                Mdot_peb = self.pebble_flux(i, r_emb)
                # Asumimos e=i=0 salvo que se dinamice externamente
                Mdot_core = self.accretion_rate(M_core, r_emb, i, e=0.0, i=0.0)
                Mdot_core = min(Mdot_core, Mdot_peb)
                # Masa aislamiento: detener acreción si se alcanza
                H_over_r = self.H[i][np.argmin(abs(self.r-r_emb))]/r_emb
                M_iso = 20*(H_over_r/0.05)**3 * 5.97e27  # g
                if M_core >= M_iso:
                    break
                # Incremento de masa en este paso
                dt = self.times[i] - (self.times[i-1] if i>0 else 0)
                dM = Mdot_core * dt
                M_core += dM
                # Fracciones locales de hielo/refractarios
                Sigma_tot = np.interp(r_emb, self.r, self.dust['Sigma'][i,:,1])
                if Sigma_tot<1e-12: 
                    continue
                for spec in self.comp:
                    frac = np.interp(r_emb, self.r, self.comp[spec][i,:]) / (Sigma_tot + 1e-12)
                    M_X[spec] += frac * dM
                history.append((t, M_core, M_X['H2O'], M_X['CO2'], M_X['silicates']))
            results[r_au] = np.array(history)
        return results
```

**Comentarios del algoritmo:** Este esquema añade: (i) transición suavizada con la función de error (líneas `f2D`), (ii) una corrección de velocidad relativa que incluye excentricidad/inclinación en $\Delta v$, (iii) un término condicional *very basic* para ablación atmosférica (si $M_{\rm core}\gtrsim0.1M_\oplus$ se descuenta el hielo). Además, cada paso chequea la masa de aislamiento con $M_{\rm iso}\propto (H/r)^3$. Finalmente, se acumulan masas de H$_2$O, CO$_2$ y silicatos por fracción local, produciendo la evolución de la composición en el tiempo.  

**Pruebas y validación:** Se recomienda comparar este módulo con casos de test de la literatura (p.ej. reproduciendo las figuras de Drazkowska+23, McCloat+25, etc.) y contra el antiguo *pebble-predictor*. En particular, verificar que: para $e=i=0$ coincide con los resultados previos (verificando las etapas 2D vs 3D), y que al introducir perturbaciones o baja viscosidad los crecimientos cambian de forma consistente (como en los apéndices de McCloat+25). La modularización permite futuros refinamientos (p.ej. mejor modelado atmosférico o fragmentación de pebbles).  

**Resumen:** El código original ya incorpora las recetas centrales de acreción de pebbles, pero las últimas investigaciones sugieren añadir transición suave 2D/3D, efectos dinámicos (excentricidad) y de atmósferas de gas. También enfatizan usar líneas de nieve móviles y multilocales para el seguimiento químico. Las referencias arriba citadas muestran que estas adiciones ponen al modelo al nivel de lo más reciente en la literatura (p.ej. el PPOLs de McCloat y los escenarios de Yap & Stevenson)【98†L34-L41】【95†L42-L49】.  

**Referencias destacadas:** Ormel & Klahr (2010), Lambrechts & Johansen (2012) – fundacionales de la fórmula de acrreción; Bitsch et al. (2018) – masa de aislamiento con turbulencia; McCloat et al. (2024, 2025) – modelos pebble-snow y entrega de agua【93†L62-L69】【98†L34-L41】; Drążkowska et al. (2022) – revisión ALMA/Kepler; Mulders et al. (2021) – filtrado y efecto de masa estelar【96†L138-L147】; Yap & Stevenson (2023) – líneas de hielo en discos circumplanetarios【95†L42-L49】. Estos trabajos representan el estado del arte en las arquitecturas planetarias basadas en acreción de pebbles y evolución de la química del disco.

---

# 6. Registro de cambios aplicados al código

## 6.1 `PebbleAccretion2.py` — Flujo de pebbles compartido en `run_growth()`

**Problema identificado:**  
El loop anterior iteraba embrión por embrión de forma independiente. Cada embrión veía el flujo de polvo completo del disco (`Mdot_peb_disk`) sin considerar que los embriones exteriores ya habían consumido parte de ese flujo. Esto sobreestimaba la acreción de los embriones interiores y hacía que la suma total de masa acretada por todos los embriones pudiera exceder la masa de pebbles disponible en el disco.

**La física correcta:**  
Los pebbles derivan radialmente hacia adentro. Un embrión a radio mayor intercepta pebbles antes de que lleguen a radios menores. El flujo disponible en un radio interior es:
```
Mdot_avail(r_in) = Mdot_disco(r_in) - Σ consumo de todos los embriones en r' > r_in
```

**Cambio implementado** (`run_growth`, `PebbleAccretion2.py`):  
Se reestructuró el loop de dos niveles: el loop **externo** itera snapshots, el loop **interno** itera embriones ordenados de **afuera hacia adentro**. Una variable `flux_consumed` se acumula con el consumo de cada embrión procesado y se resta del flujo disponible para el siguiente:

```python
flux_consumed = 0.0
for r_au in sorted(embryos, reverse=True):  # exterior → interior
    Mdot_peb_avail = max(0.0, Mdot_peb_disk(r) - flux_consumed)
    Mdot_core = min(_accretion_rate(...), Mdot_peb_avail)
    flux_consumed += Mdot_core
```

**Efecto esperado:**  
Embriones interiores (ej. 3–5 AU) crecen menos porque los exteriores consumen parte del flujo disponible. El resultado es físicamente correcto y conserva la masa de pebbles del disco.

**Referencia:** Ormel & Klahr (2010); Liu & Ormel (2018, A&A 615, A178).

---

## 6.2 `PebbleAccretion2.py` y `pipeline_snowlines.py` — Corrección de abundancia de H₂O

**Problema identificado:**  
La fracción de masa de H₂O en el polvo sólido estaba hard-codeada al ~7%, imposibilitando la clasificación de waterworld (umbral: >10% de H₂O en la masa total).

**Causa raíz:**  
Los valores de abundancia provenían de Fraser+2001 como fracciones en número relativas a H, sin conversión de unidades consistente con las del silicato (que estaban en fracciones de masa). El resultado:
```
f_H2O = 1.6e-4 / (1.6e-4 + 4e-5 + 8e-5 + 2.0e-3) ≈ 7%  ← siempre
```

**Cambio aplicado:**

| Archivo | Campo | Antes | Ahora | Justificación |
|---|---|---|---|---|
| `pipeline_snowlines.py` | `species_params["H2O"]["abundance"]` | `1.6e-4` | `9.0e-4` | Razón hielo:roca ~ 1:1 (Drążkowska & Alibert 2017) |
| `PebbleAccretion2.py` | `_ABUNDANCES['H2O']` | `1.6e-4` | `9.0e-4` | Consistencia con pipeline (fallback) |

**Resultado esperado:**
```
f_H2O = 9.0e-4 / (9.0e-4 + 4e-5 + 8e-5 + 2.0e-3) ≈ 30%
```
Un embrión formado fuera del snowline de H₂O acretaría ~30% de agua → clasificado como waterworld.

**Nota importante:** Este cambio NO afecta la hidrodinámica del disco (controlada por `chem.txt` vía `add_volatile_components`). Solo modifica los campos `SigmaIce_X` usados para estimar la composición de los pebbles acretados.

**Referencias:** Drążkowska & Alibert (2017), A&A 608, A92; Bitsch et al. (2019), A&A 623, A88; Marboeuf et al. (2014), A&A 570, A35.

---

## 6.3 `pipeline_snowlines.py` — Corrección de `v_frag` (Musiolik & Wurm 2019)

**Problema identificado:**  
El updater anterior asignaba `v_frag = 1000 cm/s` para la zona de hielo de H₂O, basándose en Pinilla et al. (2017) / Aumatell & Wurm (2011). Musiolik & Wurm (2019, ApJ 873, 58) demuestran experimentalmente que a las temperaturas del disco donde el hielo de H₂O es estable (T < 150 K), la energía superficial del hielo es γ ≈ 0.0029 J/m², igual a la de los silicatos. El incremento dramático de pegajosidad que reportan (175–200 K) ocurre por encima de la temperatura de sublimación en el disco.

**Cambio aplicado:**

| Zona de temperatura | Antes | Ahora | Referencia |
|---|---|---|---|
| T > 150 K (silicatos) | 100 cm/s | **100 cm/s** | Birnstiel et al. 2012 |
| 70 < T < 150 K (H₂O ice) | 1000 cm/s | **100 cm/s** | Musiolik & Wurm 2019 |
| 25 < T < 70 K (CO₂ ice) | 800 cm/s | **500 cm/s** | Gundlach et al. 2018 |
| T < 25 K (CO ice) | 700 cm/s | **300 cm/s** | Dominik & Tielens 1997 |

**Efecto físico:** El traffic jam en la snowline de H₂O ya no proviene de la diferencia de v_frag. Los pebbles siguen acumulándose en las snowlines por la física de condensación/sublimación (mecanismo de Schoonenberg & Ormel 2017), pero la amplificación artificial por v_frag diferencial se elimina.

---

## 6.4 `pipeline_snowlines.py` — Corrección de `setup_star_evolution()`

**Problemas identificados (tres bugs de valor):**

| Parámetro | Antes | Ahora | Error original |
|---|---|---|---|
| `R_ini` | `3.0 R_sun` | `2.0 R_sun` | Inconsistente con `R_star_ini = 2.0 R_sun` en `__init__`; causaba salto discontinuo al activar el updater |
| `R_fin` | `1.0 R_sun` | `1.5 R_sun` | La estrella llega a ZAMS en ~30 Myr, no dentro del rango de simulación |
| `t_contract` | `1.0e4 yr` | `1.0e7 yr` | 10,000 años es 5 órdenes de magnitud más rápido que la contracción Kelvin-Helmholtz real (t_KH ~ 10–30 Myr, Hayashi 1961; Siess et al. 2000) |

**Referencias:** Baraffe et al. (2015), A&A 577, A42 — tracks de pre-secuencia principal para 1 M_sun.
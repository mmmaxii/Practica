# 1. Física de la acreción de pebbles 

La acreción de pebbles (guijarros) sobre embriones planetarios se describe analíticamente diferenciando dos regímenes según la verticalización del disco de pebbles: el **2D** (pebbles concentrados cerca del plano) y el **3D** (distribuidos según una escala de altura $H_{\rm peb}$). El tamaño de captura efectivo (o radio de colisión) del embrión depende de la masa $M$ del núcleo, la viscosidad del gas y el tamaño aerodinámico de los pebbles (Stokes $St$). En general se definen dos radios característicos: el **radio de Hill** 
$$R_H=r_p\left(\frac{M}{3M_*}\right)^{1/3},$$ 
y el **radio de Bondi** 
$$R_B=\frac{2GM}{(\Delta v)^2},$$ 
donde $\Delta v\sim \eta v_K$ es la velocidad relativa entre el núcleo y el gas/pebbles (con $\eta$ el parámetro de gradiente de presión del gas)【32†L0-L2】. Para pebbles fuertemente acoplados (Stokes pequeño) el radio de captura efectivo suele ser $R_{\rm acc}\sim R_B\sqrt{St}$, mientras que para pebbles con desacoplamiento moderado se usa $R_{\rm acc}\sim R_H\,(St/0.1)^{1/3}$ (Lambrechts & Johansen 2012; Ormel & Klahr 2010). Una expresión unificada aproximada es tomar el mínimo entre estas escalas: 
$$R_{\rm acc}\approx \min\Big[R_H,\,R_B\,\sqrt{St}\Big].$$ 

Con $R_{\rm acc}$ se escribe la tasa de acreción del núcleo en 2D como: 
$$\dot M_{\rm core}^{(2D)} \approx 2\,R_{\rm acc}\,\Sigma_{\rm peb}\,\Delta v,$$ 
y en 3D como: 
$$\dot M_{\rm core}^{(3D)} \approx \pi R_{\rm acc}^2\,\rho_{\rm peb}\,\Delta v,$$ 
donde $\Sigma_{\rm peb}$ es la densidad de superficie de pebbles y $\rho_{\rm peb}\approx \Sigma_{\rm peb}/(\sqrt{2\pi}H_{\rm peb})$ su densidad volumétrica. El paso de 3D a 2D ocurre cuando el radio de captura supera la escala vertical de los pebbles: $R_{\rm acc}\gtrsim H_{\rm peb}$. En dicho límite $H_{\rm peb}=H\,\sqrt{\alpha/(\alpha+St)}$ para turbulencia $\alpha$. Se observa que la eficiencia de acreción crece con $M$ y con $St$, alcanzando un máximo alrededor de $St\sim 0.1$ (Ormel & Klahr 2010; Liu & Ormel 2018). En la práctica, se implementa el régimen efectivo usando la fórmula mayor entre $\dot M_{2D}$ y $\dot M_{3D}$, o con una transición suave (p.ej. Ida & Guillot 2016). 

**Bibliografía:** Ormel & Klahr (2010), Lambrechts & Johansen (2012, 2014), Liu & Ormel (2018) ofrecen derivaciones detalladas de estas fórmulas (ver especialmente Ormel & Klahr 2010 para la captura gravitatoria con gas).

# 2. Cálculo del flujo de pebbles Ṁ_pebble(r,t) 

El flujo de pebbles en un disco se define como la tasa de masa de sólidos (guijarros) que atraviesa una circunferencia de radio $r$. En CGS se calcula mediante: 
$$\dot M_{\rm pebble}(r,t)\;=\;2\pi\,r\,\Sigma_{\rm peb}(r,t)\,|v_r(r,t)|,$$ 
donde $\Sigma_{\rm peb}$ es la densidad de superficie de la población de pebbles (la componente grande del polvo) y $v_r$ su velocidad radial de deriva. En los datos de TripodPy, **dust.Sigma[t,r,1]** (o índice adecuado) corresponde a la superficie de la población de grano grande (pebbles), y **dust.v.rad[t,r,i]** son las velocidades de 5 tamaños muestreados. Una opción práctica es usar el tamaño más grande (o hacer promedio ponderado por masa) para $v_r$. Por ejemplo: 
- Tomar $\Sigma_{\rm peb}(r)=\texttt{dust.Sigma}[t,r,1]$ (distribución grande). 
- Tomar $v_r$ = valor absoluto de $\texttt{dust.v.rad}[t,r,k]$ para $k$ del mayor Stokes (o hacer $\Sigma_i \Sigma_i v_i/\Sigma_i$). 

Así se obtiene $\dot M_{\rm pebble}(r,t)$ en g/s. Este flujo se puede pre-computar para cada snapshot $t$ y radio $r$ útil. Alternativamente, puede calcularse localmente en el código: en cada paso de tiempo, extraer $\Sigma_{\rm peb}$ e interpolar $v_r$ hasta la posición del embrión. En resumen, se usa **dust.Sigma** (bin grande) y **dust.v.rad** para estimar el flujo de sólidos.

# 3. Seguimiento de la composición química 

La composición de la masa acumulada por el embrión se obtiene evaluando las fracciones de cada especie en el flujo de pebbles. Dados los campos químicos de salida (**components.X.dust.Sigma**), se procede así: para cada especie $X$ (H$_2$O, CO$_2$, silicato, etc.), definimos la fracción en el disco en $r,t$ como 
$$f_X(r,t) = \frac{\Sigma_{X}(r,t)}{\Sigma_{\rm peb}(r,t)},$$ 
donde $\Sigma_X$ es la superficie de polvo helado de esa especie (por ejemplo **components.H2O.dust.Sigma**). En cada paso de tiempo $\Delta t$, el núcleo gana $\Delta M_{\rm core}=\dot M_{\rm core}\,\Delta t$ de pebbles, de los cuales $\Delta M_X=f_X\,\Delta M_{\rm core}$ corresponde a la especie $X$. Llevando contabilidad acumulativa: 
$$M_{X}(t+\Delta t) = M_X(t) + f_X(r_{\rm emb},t)\,\Delta M_{\rm core},$$ 
para cada especie. Finalmente la fracción final $X$ en el núcleo es $M_X/M_{\rm core}$. Es esencial usar las fracciones instantáneas $f_X(r_{\rm emb},t)$ calculadas desde las salidas, de modo que cambian con el tiempo y con la posición relativa al embrión. En CGS, basta acumular masas en gramos. De esta forma, la composición del núcleo evoluciona con $\dot M_{\rm core}$ y las fracciones del disco.

# 4. Efecto del cruce de la línea de nieve 

Cuando la línea de nieve de, digamos, H$_2$O migra durante la evolución del disco, un embrión puede pasar de verse con pebbles secos a pebbles ricos en hielo (o viceversa). La implementación es sencilla: en cada instante se compara la posición fija $r_{\rm emb}$ con la línea de nieve $r_{\rm snow}(t)$ de cada especie. Por ejemplo, si $r_{\rm emb}$ estaba **dentro** de $r_{\rm snow,H2O}(t)$ (embrión seco) no recibe H$_2$O en los pebbles ($f_{H2O}=0$), pero cuando la línea de nieve se mueve por dentro de $r_{\rm emb}$, pasa a estar **afuera** de ella, de modo que $f_{H2O}>0$ según los campos de hielo. En la práctica, basta evaluar $f_X(r_{\rm emb},t)$ usando los datos de especies: si el embrión está dentro del radio de hielo (por debajo de la línea), el campo $\Sigma_X$ dado ya estará en cero (o muy bajo), y $f_X$ reflejará la transición. Es decir, **no** es necesario un manejo especial: las fracciones de hielo saltan de 0 a un valor positivo cuando $r_{\rm emb}$ cruza la línea. Lo importante es usar los valores de $r_{\rm snow}(t)$ (disponibles en *grid.rsnow_X*) para verificar dónde se sitúa el embrión cada instante. Así el modelo incorporará automáticamente que, a partir del cruce, la aportación de H$_2$O (u otras especies) cambia en el flujo de pebbles y por tanto en la composición acumulada.

# 5. Masa de aislamiento de pebbles 

La **masa de aislamiento de pebbles** $M_{\rm iso}$ es la masa del núcleo a la cual se forma un reto en la presión local que impide la deriva de pebbles adicionales. Lambrechts et al. (2014) dieron una fórmula aproximada (para un disco con caída de temperatura típica)  
$$M_{\rm iso}\approx 20\,M_\oplus\left(\frac{H/r}{0.05}\right)^3,$$ 
donde $(H/r)$ es el aspecto vertical del disco en el lugar del planeta【7†L17-L19】. Bitsch et al. (2018) refinaron esto incluyendo turbulencia y gradiente de presión:  
$$M_{\rm iso}\approx 25\,M_\oplus\left(\frac{H/r}{0.05}\right)^3\Big[0.34\Big(\log_{10}\frac{\alpha}{10^{-3}}+3\Big)+0.66\Big] \Big[1-\frac{\partial\ln P/\partial\ln r+2.5}{6}\Big],$$ 
donde $\alpha$ es la viscosidad turbulenta (estimada por $\nu/(c_s H)$) y $\partial\ln P/\partial\ln r$ proviene de la presión radial (por ejemplo de *gas.eta*). En CGS: $H=r\,c_s/v_K$ con $v_K=\sqrt{GM_*/r}$ (usando $M_*=1\,M_\odot$ si no se indica otra). Con *gas.nu*, *gas.cs* y *gas.eta* se calcula $\alpha$ y $H/r$ en cada radio. Al alcanzar $M_{\rm core}\ge M_{\rm iso}$, se asume que $\dot M_{\rm core}\to 0$ (ceasing pebble accretion).  

Por simplicidad, en código se suele usar la fórmula de Lambrechts-Johansen como cota base y opcionalmente multiplicadores por viscous; p.ej. $M_{\rm iso}\sim 25\,(H/r)^3\,M_\oplus$ como orden de magnitud. Recordar usar las unidades de masa terrestres convertidas a gramos para consistencia CGS.

# 6. Estructura de la clase Python *PebbleAccretionModule* 

A continuación se sugiere una estructura en Python (CGS) que lee snapshots HDF5 de TripodPy y calcula $M_{\rm core}(t)$ y composición para una lista de embriónes:

```python
class PebbleAccretionModule:
    def __init__(self, h5file, M_star=1.0):
        import h5py, numpy as np
        # Cargar archivo HDF5 de TripodPy
        f = h5py.File(h5file, 'r')
        self.grid = {key: f['grid'][key][:] for key in ['r','rsnow_H2O','rsnow_CO2','rsnow_CO']}
        self.gas = {key: f['gas'][key][:] for key in ['Sigma','T','eta','cs','nu']}
        self.dust = {
            'Sigma': f['dust/Sigma'][:],    # shape (t,r,2)
            'vr':    f['dust/v.rad'][:],    # shape (t,r,5)
            'St':    f['dust/St'][:]        # shape (t,r,5)
        }
        self.comp = {
            'H2O':   f['components/H2O/dust/Sigma'][:],    # (t,r)
            'CO2':   f['components/CO2/dust/Sigma'][:],
            'sil':   f['components/silicates/dust/Sigma'][:]
        }
        self.times = f['t'][:]   # tiempos en s
        self.rgrid = self.grid['r'][0,:]  # suponemos grid fijo o casi fijo
        self.M_star = M_star * 1.989e33     # masa estelar en g
        f.close()
        # Precalcular Omega y H
        G = 6.674e-8
        self.Omega = np.sqrt(G*self.M_star/self.rgrid**3)
        self.H = self.gas['cs'][0,:] / self.Omega  # escala de gas (t indep sup)
        # Calcular alpha si es necesario
        self.alpha = np.zeros_like(self.rgrid)
        # Evitar división por cero
        with np.errstate(divide='ignore', invalid='ignore'):
            self.alpha = self.gas['nu'][0,:] / (self.gas['cs'][0,:] * self.H)
            self.alpha[np.isnan(self.alpha)] = 1e-3

    def pebble_flux(self, t_idx, r_emb):
        # Interp. de Sigma_peb y v_rad en radio del embrión
        import numpy as np
        r = self.rgrid
        Sigma_large = self.dust['Sigma'][t_idx,:,1]  # elegir bin 1 para pebbles grandes
        vr_samples  = self.dust['vr'][t_idx,:,:]     # shape (r,5)
        # Elegir la velocidad del mayor St (último índice)
        v_r = vr_samples[:, -1]
        Sigma_p = np.interp(r_emb, r, Sigma_large)
        v_r_interp = np.interp(r_emb, r, v_r)
        Mdot_peb = 2*np.pi * r_emb * Sigma_p * abs(v_r_interp)
        return Mdot_peb  # g/s

    def accretion_rate(self, M_core, r_emb, t_idx):
        # Calcular eficiencia y Mdot_core usando fórmulas (2D/3D)
        import numpy as np
        # Extraer parámetros locales
        r = r_emb
        # Gas parameters en r_emb
        # Interpolar T, cs, eta, etc.
        T = np.interp(r, self.rgrid, self.gas['T'][t_idx])
        cs = np.interp(r, self.rgrid, self.gas['cs'][t_idx])
        Omega = np.interp(r, self.rgrid, self.Omega)
        eta = np.interp(r, self.rgrid, self.gas['eta'][t_idx])
        # Stokes y escala de pebbles (aprox)
        St = np.interp(r, self.rgrid, self.dust['St'][t_idx,:,4])  # mayor tamaño muestrado
        # Calcular radios relevantes
        G = 6.674e-8
        R_H = r * (M_core/(3*self.M_star))**(1/3)
        delta_v = eta * np.interp(r, self.rgrid, self.gas['nu'][t_idx]) # approx
        R_B = 2*G*M_core/(delta_v**2 + 1e-10)
        # Radio de acreción efectivo
        R_acc = min(R_H, R_B*np.sqrt(St))
        # Altura del disco de pebbles
        H_peb = self.H * np.sqrt(self.alpha[np.argmin(abs(self.rgrid-r))]/(self.alpha[np.argmin(abs(self.rgrid-r))]+St))
        # Elegir formula 2D o 3D
        Sigma_p = np.interp(r, self.rgrid, self.dust['Sigma'][t_idx,:,1])
        rho_p = Sigma_p/(np.sqrt(2*np.pi)*H_peb + 1e-10)
        if R_acc > H_peb:
            # Regímen 2D
            eps = 2 * R_acc * abs(delta_v)
            Mdot = eps * Sigma_p
        else:
            # Regímen 3D
            eps = np.pi * R_acc**2 * abs(delta_v) * rho_p
            Mdot = eps
        return Mdot  # g/s

    def run_embryo(self, r_emb_AU, M0=1e24):
        import numpy as np
        # Simular crecimiento de un embrión a r_emb (AU) 
        r_emb = r_emb_AU * 1.4959787e13  # pasar a cm
        M_core = M0  # masa inicial del núcleo en g (p.ej. muy pequeño)
        M_H2O = 0.0; M_CO2 = 0.0; M_sil = 0.0
        history = []
        for i,t in enumerate(self.times):
            # Calcular flujo de pebbles y acreción
            Mdot_peb = self.pebble_flux(i, r_emb)
            Mdot_core = self.accretion_rate(M_core, r_emb, i)
            # Evitar sobrepasar el flujo disponible
            Mdot_core = min(Mdot_core, Mdot_peb)
            # Parar si alcanza masa de aislamiento
            H_by_r = np.interp(r_emb, self.rgrid, self.H/self.rgrid)
            M_iso = 20 * (H_by_r/0.05)**3 * 5.97e27  # 20 M_earth en g, Lambrechts14
            if M_core >= M_iso:
                break
            # Actualizar masas
            dt = (self.times[i] - self.times[i-1]) if i>0 else self.times[0]
            dM = Mdot_core * dt
            # Fracciones de especies locales
            Sigma_tot = np.interp(r_emb, self.rgrid, self.dust['Sigma'][i,:,1])
            # Si Sigma_tot ~0, salteo
            if Sigma_tot < 1e-30:
                continue
            f_H2O = np.interp(r_emb, self.rgrid, self.comp['H2O'][i,:]) / Sigma_tot
            f_CO2 = np.interp(r_emb, self.rgrid, self.comp['CO2'][i,:]) / Sigma_tot
            f_sil = np.interp(r_emb, self.rgrid, self.comp['sil'][i,:]) / Sigma_tot
            # Acumular
            M_H2O += f_H2O * dM
            M_CO2 += f_CO2 * dM
            M_sil += f_sil * dM
            M_core += dM
            # Guardar histórico
            history.append((t, M_core, M_H2O, M_CO2, M_sil))
        return np.array(history)
```

Este esquema lee los datos de TripodPy, calcula paso a paso la acreción y va acumulando $M_{\rm core}$ junto con la masa de H$_2$O, CO$_2$ y silicato. Se pueden adaptar las fórmulas de $\dot M_{\rm core}$ en `accretion_rate()` para incluir con detalle los regímenes 2D/3D exactos (p.ej. usando la dependencia completa con Stokes o funciones de eficiencia). En cualquier caso, el módulo propuesto es compatible con los outputs HDF5 de TripodPy, usa unidades CGS, y produce como salida la evolución temporal de la masa del núcleo y su composición química acumulada.

**Fuentes:** Las fórmulas e ideas de acreción de pebbles siguen modelos estándar (Ormel & Klahr 2010; Lambrechts & Johansen 2012, 2014; Liu & Ormel 2018) y los conceptos de masa de aislamiento (Lambrechts et al. 2014; Bitsch et al. 2018). Las expresiones numéricas y el algoritmo fueron adaptados al formato de salida de TripodPy mencionado. Cada fórmula clave está tomada de estas referencias clásicas en el estudio de acreción de pebbles.
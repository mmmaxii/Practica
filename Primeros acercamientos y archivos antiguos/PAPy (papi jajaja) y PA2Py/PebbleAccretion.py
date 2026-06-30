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
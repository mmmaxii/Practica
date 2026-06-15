import os
import glob
import pickle
import numpy as np
import sys
import contextlib

# Asegurarse de poder importar PA3Py
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
try:
    from PA3Py.PebbleAccretion3 import PebbleAccretionModule3
    import dustpy.constants as c
except ImportError as e:
    print(f"Advertencia: No se pudo importar PA3Py/dustpy ({e}).")

class BaseAnalyzer:
    def __init__(self, runs_dir, cache_file, m0_earth, alpha_val=None):
        self.runs_dir = runs_dir
        self.cache_file = cache_file
        self.m0_earth = m0_earth
        self.alpha_val = alpha_val
        self.data = None

    def get_target_runs(self):
        runs = glob.glob(os.path.join(self.runs_dir, "*"))
        target_runs = []
        for rpath in runs:
            if not os.path.isdir(rpath):
                continue
                
            run_name = os.path.basename(rpath)
            if self.alpha_val is not None:
                # Asegurarse de que coincida el alpha exacto (evitar _a0.001 en _a0.0015)
                if not (run_name.endswith(f"_a{self.alpha_val}") or f"_a{self.alpha_val}_" in run_name):
                    continue
                    
            if os.path.exists(os.path.join(rpath, "data0000.hdf5")):
                completed = os.path.exists(os.path.join(rpath, "data0099.hdf5"))
                target_runs.append((rpath, completed))
        target_runs.sort()
        return target_runs

    def parse_run_name(self, run_name):
        raise NotImplementedError("Las subclases deben implementar parse_run_name")

    def extract_data(self):
        runs_info = self.get_target_runs()
        self.data = []
        
        print(f"\nIniciando extracción PA3 para {len(runs_info)} simulaciones en {self.runs_dir}...")
        for i, (rpath, completed) in enumerate(runs_info):
            run_name = os.path.basename(rpath)
            try:
                parsed_params = self.parse_run_name(run_name)
            except Exception as e:
                # Si falla el parseo, probablemente no es una carpeta de run válida
                continue
                
            print(f"[{i+1}/{len(runs_info)}] Procesando {run_name} (Completado: {completed})...")
            try:
                with open(os.devnull, 'w') as fnull, contextlib.redirect_stdout(fnull):
                    pa3 = PebbleAccretionModule3.from_datadir(rpath, M_star=1.0)
                    res = pa3.run_growth([1.0], M0_g=self.m0_earth * c.M_earth)
            except Exception as e:
                print(f"  Error cargando HDF5 en {run_name}: {e}")
                continue
            hist = res[1.0]
            
            if len(hist) == 0:
                print(f"  Sin evolución registrada para {run_name}")
                continue
                
            times_yr = hist[:, 0] / c.year
            mass_e = hist[:, 1] / c.M_earth
            m_h2o = hist[:, 2]
            m_sil = hist[:, 3]
            rsnow_au = hist[:, 4]
            m_iso_e = hist[:, 5] / c.M_earth
            
            total_final = m_h2o[-1] + m_sil[-1]
            frac_h2o_final = 100.0 * (m_h2o[-1] / total_final) if total_final > 0 else 0.0
            
            cross_idx = np.where(rsnow_au < 1.0)[0]
            t_cross_1au = times_yr[cross_idx[0]] if len(cross_idx) > 0 else None
            
            entry = {
                'run_name': run_name,
                'times_yr': times_yr,
                'mass_e': mass_e,
                'm_iso_e': m_iso_e,
                'frac_h2o_final': frac_h2o_final,
                't_cross_1au': t_cross_1au,
                'completed': completed,
                'm0_earth': self.m0_earth
            }
            entry.update(parsed_params)
            self.data.append(entry)
            
        os.makedirs(os.path.dirname(self.cache_file), exist_ok=True)
        with open(self.cache_file, 'wb') as f:
            pickle.dump(self.data, f)
            
        return self.data

    def load_or_extract(self):
        if os.path.exists(self.cache_file):
            print(f"Cargando datos desde caché: {self.cache_file}")
            with open(self.cache_file, 'rb') as f:
                self.data = pickle.load(f)
        else:
            self.extract_data()
        return self.data


class GeneralAnalyzer(BaseAnalyzer):
    def parse_run_name(self, run_name):
        # Ej: run_r10.0_m0.05_a0.001 o run_r10.0_m0.05_a0.001_v1.0
        parts = run_name.split('_')
        r_val = float(parts[1][1:])
        m_val = float(parts[2][1:])
        a_val = float(parts[3][1:])
        return {'r_gap': r_val, 'M_gap': m_val, 'alpha': a_val}


class SinusoidalAnalyzer(BaseAnalyzer):
    def parse_run_name(self, run_name):
        parts = run_name.split('_')
        if len(parts) > 1 and "ngap" in parts[1]:
            # Nuevo formato Geryon: run_ngap5_A1.0_a0.0005_rmin0.7
            n_gaps = int(parts[1].replace('ngap', ''))
            amp_val = float(parts[2][1:])
            alpha_val = float(parts[3][1:])
            
            if n_gaps <= 3:
                gap_type = "PocosGaps"
            elif n_gaps <= 5:
                gap_type = "Intermedio"
            else:
                gap_type = "MuchosGaps"
                
            return {'gap_type': gap_type, 'n_gaps': n_gaps, 'amp_val': amp_val, 'alpha': alpha_val}
        else:
            # Formato viejo local: PocosGaps_A0.5
            gap_type = parts[0]
            amp_val = float(parts[1][1:])
            return {'gap_type': gap_type, 'amp_val': amp_val, 'alpha': 0.001}


class DelayedAnalyzer(BaseAnalyzer):
    def parse_run_name(self, run_name):
        # Ej: run_retrasado_r10.0_m0.01_tgap1000000
        parts = run_name.split('_')
        r_val = float(parts[2][1:])
        m_val = float(parts[3][1:])
        tgap_str = parts[4].replace('tgap', '')
        tgap_val = float(tgap_str)
        return {'r_gap': r_val, 'M_gap': m_val, 't_gap': tgap_val}

class SmoothAnalyzer(BaseAnalyzer):
    def parse_run_name(self, run_name):
        # Ej: run_smooth_a0.0001_v10
        parts = run_name.split('_')
        a_val = float(parts[2][1:])
        return {'alpha': a_val}


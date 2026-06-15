import os
import glob
import h5py
import numpy as np

class SnapshotAnalyzer:
    """
    Analizador de alto rendimiento para extraer la hidrodinámica 
    espacio-temporal a partir de snapshots HDF5 de TripodPy.
    """
    def __init__(self, run_dir):
        self.run_dir = run_dir
        self.h5_files = sorted(glob.glob(os.path.join(run_dir, "*.hdf5")))
        if not self.h5_files:
            raise FileNotFoundError(f"No se encontraron archivos HDF5 en {run_dir}")
            
        # Extraer parámetros constantes del grid usando el primer snapshot
        with h5py.File(self.h5_files[0], 'r') as f:
            self.r = f['grid/r'][:]   # Centros de celda [cm]
            self.ri = f['grid/ri'][:] # Interfaces de celda [cm]
            
        self.AU = 1.495978707e13 # 1 AU en cm
        self.idx_1au_r = np.argmin(np.abs(self.r - self.AU))
        self.idx_1au_ri = np.argmin(np.abs(self.ri - self.AU))

    def extract_spatiotemporal_data(self, subsampling=1):
        """
        Extrae el tiempo, la relación polvo-gas (epsilon), el flujo local de pebbles a 1 AU
        y el tamaño máximo de grano (amax) a 1 AU.
        """
        files_to_process = self.h5_files[::subsampling]
        N_times = len(files_to_process)
        N_r = len(self.r)
        
        t_array = np.zeros(N_times)
        eps_matrix = np.zeros((N_times, N_r))
        flux_1au = np.zeros(N_times)
        amax_1au = np.zeros(N_times)
        
        for idx, file_path in enumerate(files_to_process):
            with h5py.File(file_path, 'r') as f:
                t_array[idx] = f['t'][()] # Tiempo actual en segundos
                
                # Perfiles de densidad superficial
                sigma_g = f['gas/Sigma'][:]
                # Sumar bines de polvo (small + large en TriPoD)
                sigma_d = np.sum(f['dust/Sigma'][:], axis=-1) 
                
                # Relación Polvo-Gas (evitando división por cero)
                eps_matrix[idx, :] = np.where(sigma_g > 0, sigma_d / sigma_g, 0.0)
                
                # Flujo neto de pebbles cruzando la interfaz de 1 AU [g/s]
                flux_1au[idx] = 2 * np.pi * self.ri[self.idx_1au_ri] * np.sum(f['dust/Fi/tot'][self.idx_1au_ri, :])
                
                # Tamaño máximo de grano en 1 AU [cm]
                # dust/a tiene forma (N_r, 5). El último elemento es a_max.
                amax_1au[idx] = f['dust/a'][self.idx_1au_r, -1]
                
        return t_array, self.r / self.AU, eps_matrix, flux_1au, amax_1au

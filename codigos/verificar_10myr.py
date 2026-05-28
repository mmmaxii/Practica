# -*- coding: utf-8 -*-
import os
import glob
import h5py
M_sun = 1.9884e33
R_sun = 6.957e10
au = 1.496e13

base_dir = "data/runs/10Myr"
runs = glob.glob(os.path.join(base_dir, "run_*"))
runs.sort()

print("="*80)
print("REPORTE DE VERIFICACIÓN 10 Myr")
print("="*80)

completados = 0
totales = len(runs)

for run_path in runs:
    run_name = os.path.basename(run_path)
    
    # Check completeness
    file_99 = os.path.join(run_path, "data0099.hdf5")
    completado = os.path.exists(file_99)
    
    # Check physics
    file_00 = os.path.join(run_path, "data0000.hdf5")
    if not os.path.exists(file_00):
        print(f"[{run_name}] ERROR: No se encontró data0000.hdf5")
        continue
        
    try:
        with h5py.File(file_00, "r") as f:
            M = f['star']['M'][()] / M_sun
            R = f['star']['R'][()] / R_sun
            T = f['star']['T'][()]
            r_in = f['grid']['r'][0] / au
            r_out = f['grid']['r'][-1] / au
            nr = len(f['grid']['r'])
            
            # Tolerancias (r_in se lee como ~0.504 por ser el centro de la celda)
            m_ok = abs(M - 1.0) < 0.05
            r_ok = abs(R - 2.1) < 0.05
            t_ok = abs(T - 4000.0) < 1.0
            r_in_ok = abs(r_in - 0.5) < 0.1
            r_out_ok = abs(r_out - 100.0) < 5.0
            nr_ok = abs(nr - 300) < 50  # 300 base + celdas refinadas
            
            if m_ok and r_ok and t_ok and r_in_ok and r_out_ok and nr_ok:
                phys_status = "Física OK"
            else:
                phys_status = f"Física ANÓMALA (M={M:.2f}, R={R:.2f}, T={T:.0f}, r=[{r_in:.2f},{r_out:.2f}], Nr={nr})"
    except Exception as e:
        phys_status = f"Error leyendo HDF5: {e}"
        
    if completado:
        status = "COMPLETADO (10 Myr)"
        completados += 1
    else:
        status = "INCOMPLETO"
        
    print(f"[{run_name:30s}] | {status:20s} | {phys_status}")
    
print("-" * 80)
print(f"Total completados: {completados}/{totales}")
print("="*80)

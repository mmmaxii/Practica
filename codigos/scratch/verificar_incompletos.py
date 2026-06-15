# -*- coding: utf-8 -*-
import os
import glob

def auditar_runs(base_dir="data/runs"):
    print(f"Auditando directorio: {base_dir} ...\n")
    
    # Buscar todas las carpetas que contengan al menos data0000.hdf5
    # Esto indica que es una simulación de TripodPy
    todas_las_runs = glob.glob(os.path.join(base_dir, "**", "data0000.hdf5"), recursive=True)
    
    incompletas = []
    totales = len(todas_las_runs)
    
    for run_file in todas_las_runs:
        run_dir = os.path.dirname(run_file)
        
        # Contar archivos HDF5
        hdf5_files = glob.glob(os.path.join(run_dir, "data*.hdf5"))
        num_snapshots = len(hdf5_files)
        
        # Criterio: menos de 97 snapshots (0 a 96)
        if num_snapshots < 97:
            # Obtener el nombre base del run
            run_name = os.path.basename(run_dir)
            # Extraer subcarpeta (ej: vf_10ms/general)
            rel_path = os.path.relpath(run_dir, base_dir)
            incompletas.append({
                'path': rel_path,
                'name': run_name,
                'snaps': num_snapshots
            })
            
    print(f"Total de simulaciones escaneadas: {totales}")
    print(f"Simulaciones incompletas (< 97 snapshots): {len(incompletas)}\n")
    
    if incompletas:
        print("=== LISTA DE SIMULACIONES INCOMPLETAS ===")
        # Ordenar por cantidad de snapshots (los que murieron más temprano primero)
        incompletas.sort(key=lambda x: x['snaps'])
        for inc in incompletas:
            print(f"Snapshots: {inc['snaps']:02d} | Run: {inc['path']}")
    else:
        print("¡Todas las simulaciones llegaron a 97 o más snapshots!")

if __name__ == "__main__":
    auditar_runs()

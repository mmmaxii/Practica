"""
inspect_hdf5.py
Muestra las claves internas de un archivo HDF5 de TripodPy.
"""
import h5py

path = r"C:\astro\Codigos practica + docs + papers\codigos\data\runs\10Myr\run_r10.0_m0.1_a0.001\data0099.hdf5"

def print_tree(name, obj):
    indent = "  " * name.count("/")
    if isinstance(obj, h5py.Dataset):
        print(f"{indent}{name}  -> shape={obj.shape}, dtype={obj.dtype}")
    else:
        print(f"{indent}{name}/")

with h5py.File(path, "r") as f:
    print("=== Claves raíz ===")
    f.visititems(print_tree)

import h5py
import sys

file_path = "data/runs/10Myr/run_r10.0_m1.0_a0.001/data0000.hdf5"
try:
    with h5py.File(file_path, "r") as f:
        print("Star M:", f['star']['M'][()])
        print("Star R:", f['star']['R'][()])
        print("Star T:", f['star']['T'][()])
        print("Grid r[0]:", f['grid']['r'][0])
        print("Grid r[-1]:", f['grid']['r'][-1])
        print("Grid Nr:", len(f['grid']['r']))
except Exception as e:
    print(f"Error: {e}")

import pickle
import numpy as np

with open('data/runs/10Myr_pa3_cache.pkl', 'rb') as f:
    data = pickle.load(f)

for item in data.get(0.001, []):
    if item['M_gap'] == 0.01:
        print(f"Run: {item['run_name']}")
        print(f"  Mass (first 5): {item['mass_e'][:5]}")
        print(f"  Mass (last 5): {item['mass_e'][-5:]}")
        print(f"  Length: {len(item['mass_e'])}")
        t_myr_log = np.log10(item['times_yr'] / 1e6)
        mass_e_log = np.log10(item['mass_e'])
        mask = np.isfinite(t_myr_log) & np.isfinite(mass_e_log)
        print(f"  Valid mask length: {np.sum(mask)}")

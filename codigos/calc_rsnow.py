import numpy as np

# from oka_interpolation
def mdot_time(t_myr, eta=1.5):
    mdot_1myr = 1e-8
    t_myr = np.maximum(t_myr, 1e-6) 
    return mdot_1myr * (t_myr)**(-eta)

# Mocking interpolation
mdot = 1e-7
r_snow_1e7 = 3.8
# r = 3.8 * (Mdot / 1e-7)**0.41

print("Tiempo [yr]  | Mdot [Msun/yr] | r_snow [AU] aprox")
for t_yr in [100, 215, 464, 1000, 2154, 4641, 10000, 21544, 46415, 100000]:
    t_myr = t_yr / 1e6
    m = mdot_time(t_myr)
    r = 3.8 * (m / 1e-7)**0.41
    print(f"{t_yr:10.2f} | {m:14.2e} | {r:10.2f}")

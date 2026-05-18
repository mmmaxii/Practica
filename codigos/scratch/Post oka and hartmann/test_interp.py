import numpy as np
from scipy.interpolate import interp1d

x = np.array([-12, -11, -10, -9, -8, -7])
y = np.array([0, 1, 2, 3, 4, 5])
f = interp1d(x, y, kind='linear', bounds_error=False, fill_value="extrapolate")
print(f([-6, -5, -4]))

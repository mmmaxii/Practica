from tripodpy import Simulation
import dustpy.constants as c
import numpy as np
import os

sim = Simulation()
sim.ini.grid.Nr = 150
sim.makegrids()
sim.initialize()

sigma_base = sim.gas.Sigma.copy()

def forzar_sigma_bump(sim):
    A1 = 5.0
    r0_1 = 30 * c.au  ;  width_1 = 5 * c.au
    r = sim.grid.r
    bump1 = A1 * np.exp(- (r - r0_1)**2 / (2 * width_1**2))
    sim.gas.Sigma[:] = sigma_base * (1 + bump1)
    return 0.

sim.gas.addfield("imposicion_densidad", 0.)
sim.gas.imposicion_densidad.updater.updater = forzar_sigma_bump
sim.gas.updateorder.insert(0, "imposicion_densidad")

sim.gas.imposicion_densidad.update()
sim.update()
sim.t.snapshots = np.linspace(0, 1e4*c.year, 3)
sim.writer.datadir = "bumps/pb_v7_test"
sim.writer.overwrite = True
sim.run()

from dustpy import hdf5writer
wrtr = hdf5writer()
wrtr.datadir = "bumps/pb_v7_test"
data = wrtr.read.all()
idx = np.argmin(np.abs(sim.grid.r - 30*c.au))
print("Sigma initial at 30 AU:", data.gas.Sigma[0, idx])
print("Sigma final at 30 AU:", data.gas.Sigma[-1, idx])

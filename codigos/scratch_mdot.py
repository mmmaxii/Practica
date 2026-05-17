import sys
import numpy as np
import dustpy.constants as c
from tripodpy import Simulation

sim = Simulation()
sim.initialize()

print("Gas fields:", dir(sim.gas))
print("Is Mdot there?", 'Mdot' in dir(sim.gas))
print("Is v there?", 'v' in dir(sim.gas))
if 'v' in dir(sim.gas):
    print("v fields:", dir(sim.gas.v))

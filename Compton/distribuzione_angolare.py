import comptonFunction
from setup import Detector
import numpy as np
import matplotlib.pyplot as plt
import random

myScintillator = Detector()
melec = 511
e_cs = 662
nevents = 16000
angle = np.zeros(nevents)
i = 0


while i < nevents:

    pin = np.array([e_cs, e_cs, 0, 0])
    s = myScintillator.scattering(0.5)
    p1, pout = myScintillator.acquirePhoton(pin, s)

    trigger = myScintillator.check_momentum(pout) and myScintillator.check_energy(p1, pout)
    if trigger:
        angle[i] = myScintillator.angularSampling(pout)
        i += 1

plt.figure()
ydata, edges, _ = plt.hist(angle, bins=50, label='Angular distribution')
xdata = 0.5 * (edges[1:] + edges[:-1])
plt.errorbar(xdata, ydata, np.sqrt(ydata), fmt='.', color='black', label ='Data (Poissonian Error')
plt.legend(loc='best')
plt.xlabel('Angle [°]')
plt.ylabel('Counts')
plt.show()
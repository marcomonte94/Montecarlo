import numpy as np
from particles import Particle, Material
from LinearEnergyTransfer import bethe_bloch
from braggCurve import *
from matplotlib import pyplot as plt

water = Material('Water', 8, 18, 1, 75e-6)
dx = 0.005 #cm
npart = 200 # n° of simulated particles
r = []
np.random.seed(80)

''' Same workflow as before (I've repeated myself, I know...)
'''

for i in range (npart):
    print(npart - i)

    p = Particle('Proton', 938., 1, 50)
    ein = p.eKin

    # Beam dispersion parameters (default values: 0.5, 0.001, 0.001)  
    pSigma, xSigma, ySigma = 0,0,0

    p.eKin = np.random.normal(p.eKin, pSigma)

    
    x, y = np.random.normal(0, xSigma), np.random.normal(0, ySigma)
    x, y = 0, 0
    
    xpath, ypath = [], []
    eLoss = []

    while p.eKin > 0:

        let = bethe_bloch(p, water)
        theta = samplingTheta(p, water, dx)
        de = enStraggling(p, water, let, dx)
        p.eKin = p.eKin - de
        x = x + dx* np.cos(theta)
        #y = y + dx* np.sin(theta)

        eLoss.append(let)
        xpath.append(x)
        #ypath.append(y)

    r.append(x)

r = np.array(r)

print('Energy: {} MeV'.format(ein))
print('Bragg peak position: {}'.format(r.mean()))
print('Bragg peak std: {}'.format(r.std()))
plt.figure()
plt.title('Range distribution: n°protons = 900, E = 100 MeV')
plt.hist(r, bins=30)
plt.xlabel('Range [cm]')
plt.ylabel('Events [a.u.]')
plt.show()

from particles import Particle, Material
import numpy as np
from matplotlib import pyplot as plt

p = Particle('Positron', 511, +1)
p.eKin = np.linspace(0.1, 635, 1000) 


def beta_spectrum():

    def fermi_function(p, Z_daughter):
        alpha = 1 / 137
        n = -Z_daughter * alpha * p.energy() / (p.momentum())
        return 2 * np.pi * n / (1 - np.exp(-2 * np.pi * n))

    Emax = 635 
    w = 1 + p.eKin / p.mass
    w0 = 1 + Emax / p.mass
    momentum = np.sqrt(w**2 - 1)
    spectra = fermi_function(p, 8) * (p.momentum()/p.mass) * w * ((w0 - w)**2)
    
    return spectra
#print(Emax - p.gamma())

plt.plot(p.eKin, spectra)
plt.show()
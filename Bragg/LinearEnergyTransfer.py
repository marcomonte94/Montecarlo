from particles import Particle, Material
import numpy as np
from matplotlib import pyplot as plt

''' Energy in Mev, distance in g / cm^2
'''
C = 0.307 # MeV cm ^2 / g
melec = 0.511

def bohr(p, medium):

    C = 0.307 # MeV cm ^2 / g

    let = (C * medium.Z/medium.A * medium.density * p.charge / p.beta()**2) * (np.log(2*melec*(p.beta()*p.gamma()**2)/medium.I))

    return let

def bethe_bloch(p, medium):


    C = 0.307 # MeV cm ^2 / g

    Tmax = (2*melec*(p.beta()*p.gamma())**2)/ (1 + 2*p.gamma()*(melec/p.mass)+(melec/p.mass)**2)
    let = (C * medium.Z/medium.A * medium.density * p.charge**2 / p.beta()**2) * (0.5 * np.log(2*melec*(p.beta()*p.gamma()**2)*Tmax/(medium.I**2)) - p.beta()**2)

    return let


if __name__ == '__main__':

    ''' Define some particles (heavy charge particles)
    '''
    p = Particle('Proton', 938., 1)
    p.eKin = np.linspace(0.1, 100000, 100000)


    alpha = Particle('Alpha', 4*938, 2)
    alpha.eKin = np.linspace(0.1, 1e7, 100000)


    c = Particle('c nucleus', 12*938, 6)
    c.eKin = np.linspace(0.1, 1e7, 100000)

    o = Particle('Oxygen nucleus', 16*938, 8)
    o.eKin = np.linspace(0.1, 5e7, 100000)

    ''' Define some material
    '''
    water = Material('Water', 10, 18, 1, 75e-6)
    helium = Material('Helium', 2, 4, 1.663e-04, 41.8e-6)
    iron = Material('Iron', 26, 56, 7.874, 286e-6)
    lead = Material('Lead', 82, 207, 11.35, 823e-6)


    plt.figure(1)
    plt.title('Bohr vs Bethe-Bloch; protons in $H_2O$')
    plt.plot(p.beta()*p.gamma(),bethe_bloch(p, water), label='Bethe-Bloch')
    plt.plot(p.beta()*p.gamma(), bohr(p, water), label='Bohr')
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 100)
    plt.ylim(1.5, 1e2)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')


    plt.figure(2)
    plt.title('Different particles in $H_2O$')
    plt.plot(p.beta()*p.gamma(),bethe_bloch(p, water), label=p.name)
    plt.plot(alpha.beta()*alpha.gamma(),bethe_bloch(alpha, water), label=alpha.name)
    plt.plot(c.beta()*c.gamma(),bethe_bloch(c, water), label=c.name)
    plt.plot(o.beta()*o.gamma(), bethe_bloch(o, water), label=o.name)
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 100)
    plt.ylim(1.5, 10000.)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')


    plt.figure(3)
    plt.title('Protons in different materials')
    plt.plot(p.beta()*p.gamma(), bethe_bloch(p, water), label=water.name)
    #plt.plot(p.beta()*p.gamma(), bethe_bloch(p, helium), label=helium.name)
    plt.plot(p.beta()*p.gamma(), bethe_bloch(p, iron), label=iron.name)
    plt.plot(p.beta()*p.gamma(), bethe_bloch(p, lead), label=lead.name)
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 100)
    plt.ylim(1e-1, 1e4)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')

    plt.show()



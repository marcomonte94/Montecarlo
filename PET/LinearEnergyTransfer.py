from particles import Particle 
import numpy as np 
from matplotlib import pyplot as plt 

''' Energy in Mev, distance in g / cm^2
'''
C = 0.307 # MeY cm ^2 / g
melec = 0.511

def bohr(p, rho, I):

    beta = np.linspace(0, 0.99999999, 10000)
    let = (C * 0.5 * rho * p.charge / beta**2) * (np.log(2*melec*(beta*p.gamma(beta)**2)/I))
    plt.plot(beta*p.gamma(beta), let, label='Bohr')

def bethe_bloch(p, rho, I, label):

    beta = np.linspace(0, 0.99999999, 10000)
    Tmax = (2*melec*(beta*p.gamma(beta)**2))/ (1 + 2*p.gamma(beta)*(melec/p.mass)+(melec/p.mass)**2)
    let = (C * 0.5 * rho * p.charge / beta**2) * (0.5 * np.log(2*melec*(beta*p.gamma(beta)**2)*Tmax/I**2) - beta**2)
    
    plt.plot(beta*p.gamma(beta), let, label=label)
 

''' Define some particles (heavy charge particles)
'''
p = Particle('Proton', 938., 1)
alpha = Particle('Alpha', 4*938, 2)
carbon = Particle('Carbon nucleus', 12*938, 6)
oxigen = Particle('Oxygen nucleus', 16*938, 8)

''' Define some material densities (g / cm^3 ) and ionization energy (MeV)
'''
rho_water, I_water = 1., 12.621e-6
rho_helium, I_helium = 1.663e-04, 41.8e-6
rho_iron, I_iron = 7.874, 286e-6
rho_lead, I_lead =  11.35, 823e-6


plt.figure(1)
plt.title('Bohr vs Bethe-Bloch; protons in $H_2O$')
bethe_bloch(p, rho_water, I_water, 'Bethe-Bloch')
bohr(p, rho_water, I_water)
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, 10000)
plt.ylim(1.5, 100.)
plt.legend(loc='best')
plt.xlabel('p / M c')
plt.ylabel('Energy loss [MeV cm^2 / g]')

plt.figure(2)
plt.title('Different particles in $H_2O$')
bethe_bloch(p, rho_water, I_water, p.name)
bethe_bloch(alpha, rho_water, I_water, alpha.name)
bethe_bloch(carbon, rho_water, I_water, carbon.name)
bethe_bloch(oxigen, rho_water, I_water, oxigen.name)
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, 10000)
plt.ylim(1.5, 1000.)
plt.legend(loc='best')
plt.xlabel('p / M c')
plt.ylabel('Energy loss [MeV cm^2 / g]')

plt.figure(3)
plt.title('Prhotons in different materials')
bethe_bloch(p, rho_water, I_water, "Water")
bethe_bloch(p, rho_helium, I_helium, "Helium")
bethe_bloch(p, rho_iron, I_iron, "Iron")
bethe_bloch(p, rho_lead, I_lead, "Lead")
plt.grid()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0.1, 10000)
plt.ylim(0, 1e5)
plt.legend(loc='best')
plt.xlabel('p / M c')
plt.ylabel('Energy loss [MeV cm^2 / g]')

plt.show()



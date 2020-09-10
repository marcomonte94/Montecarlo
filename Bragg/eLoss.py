from particles import Particle, Material
import numpy as np
from matplotlib import pyplot as plt

melec = 0.511

def bethe_bloch(p, medium):


    C = 0.1536 # MeV cm ^2 / g
    tau = p.eKin / melec

    def F(t):

        if p.charge < 0:
            print('Sono un elettrone')
            return 1 - p.beta()**2 + ((t**2/8) - (2*t + 1)*np.log(2)) / (t + 1)**2
        else:
            print('Sono un positrone')
            return 2*np.log(2) - ((p.beta()**2)/12) * (23 + 14/(t+2) + 10/((t+2)**2) + 4/((t+2)**3))


    let = (C * medium.Z/medium.A * medium.density * p.charge**2 / p.beta()**2) * (np.log(  ((tau+2)*tau**2) / (2*(medium.I/melec)**2)  ) + F(tau))

    return let

def bremsstrahlung(p, medium):

    re = 2.82e-13
    Na = 6.022e23
    alpha = 1 / 137

    N = medium.density * Na / medium.A

    let = []
    
    for i in p.energy():

        if i > p.mass:
            phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(2*i/p.mass) - 1/3)
        else:
            phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(183*medium.Z**(-1/3)) + 1/18)

        let.append(N * i * phi_r)
    '''
    if p.energy() > p.mass:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(2*p.energy()/p.mass) - 1/3)
    else:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(183*medium.Z**(-1/3)) + 1/18)
    
    let = N * p.energy() * phi_r
    '''
    return let


if __name__ == '__main__':

    water = Material('Water', 10, 18, 1, 75e-6)

    e = Particle('Electron', 0.511, -1)
    e.eKin=np.linspace(0000.1, 5e3, 100000)

    p = Particle('Positron', 0.511, 1)
    p.eKin=np.linspace(0000.1, 5e3, 100000)

    plt.figure(1)
    plt.title('Collision energy loss; e$^+$ e$^-$ in $H_2O$')
    plt.plot(e.eKin, bethe_bloch(e, water), label='e$^-$')
    plt.plot(p.eKin, bethe_bloch(p, water), label='e$^+$')
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(0.1, 100)
    #plt.ylim(1.5, 10)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')
    plt.ylim(1e-1, 1e4)
    
    plt.figure(2)
    plt.title('Collision energy loss; e$^+$ e$^-$ in $H_2O$')
    plt.plot(e.eKin, bethe_bloch(e, water), label='Collisional')
    plt.plot(p.eKin, bremsstrahlung(e, water), label='Radiative')
    plt.plot(p.eKin, bethe_bloch(e, water)+bremsstrahlung(e, water))
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(0.1, 100)
    #plt.ylim(1.5, 10)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')
    

    plt.show()

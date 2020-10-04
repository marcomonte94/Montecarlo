from particles import Particle, Material
import numpy as np
from matplotlib import pyplot as plt

melec = 511

def bethe_bloch(p, medium):


    #C = 0.1536 # MeV cm ^2 / g
    C = 153.6 # keV cm ^2 / g
    tau = p.eKin / p.mass

    def F(t):

        if p.charge < 0:
            #print('Sono un elettrone')
            return 1 - p.beta()**2 + ((t**2/8) - (2*t + 1)*np.log(2)) / (t + 1)**2
        else:
            #print('Sono un positrone')
            return 2*np.log(2) - ((p.beta()**2)/12) * (23 + 14/(t+2) + 10/((t+2)**2) + 4/((t+2)**3))


    let = (C * medium.Z/medium.A * medium.density * p.charge**2 / p.beta()**2) * (np.log(  ((tau+2)*tau**2) / (2*(medium.I/melec)**2)  ) + F(tau))

    return let

def bremsstrahlung(p, medium):

    re = 2.82e-13
    Na = 6.022e23
    alpha = 1 / 137

    N = medium.density * Na / medium.A

    let = []

    #for i in range(len(p.energy())):

    if float(p.energy()) > p.mass:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(2*p.energy()/p.mass) - 1/3)
    else:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(183*medium.Z**(-1/3)) + 1/18)

    let = (N * p.energy() * phi_r)
    '''
    if p.energy() > p.mass:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(2*p.energy()/p.mass) - 1/3)
    else:
        phi_r = 4 * alpha* (medium.Z**2) * re**2 * ( np.log(183*medium.Z**(-1/3)) + 1/18)

    let = N * p.energy() * phi_r
    '''
    return let


if __name__ == '__main__':

    water = Material('Water', 7.2, 18, 1, 75e-6)

    e = Particle('Electron', 511, -1)
    e.eKin=np.linspace(100, 1e7, 100000)

    p = Particle('Positron', 511, 1)
    p.eKin=np.linspace(100, 1e7, 100000)

    print(min(bethe_bloch(p, water)))

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.title('Collision energy loss; e$^+$ e$^-$ in $H_2O$')
    plt.plot(e.eKin, bethe_bloch(e, water),'--', color='blue',  label='e$^-$')
    plt.plot(p.eKin, bethe_bloch(p, water), color='red', label='e$^+$')
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e2, 1e7)
    #plt.ylim(1.5, 10)
    plt.legend(loc='best')
    plt.xlabel('$E_{kin}$ [keV]')
    plt.ylabel('dE/dx [keV/cm]')
    plt.ylim(1e3, 2e4)

    let_br = []
    for i in range(len(e.eKin)):
        #p = Particle('Positron', 511, 1)
        p.eKin = e.eKin[i]
        let_br.append(bremsstrahlung(p, water))

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.title('Collision vs radiative energy loss; e$^+$ in $H_2O$')
    p.eKin=np.linspace(100, 1e7, 100000)
    plt.plot(p.eKin, bethe_bloch(p, water),'--', color='black', label='Collision')
    plt.plot(p.eKin, let_br, '--', color='blue', label='Radiative')
    plt.plot(p.eKin, bethe_bloch(p, water)+np.asarray(let_br), color='red', label='Total')
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(100, 1e7)
    plt.ylim(1e2, 1e6)
    plt.xlabel('$E_{kin}$ [keV]')
    plt.ylabel('dE/dx [keV/cm]')
    plt.legend(loc='best')

    plt.show()

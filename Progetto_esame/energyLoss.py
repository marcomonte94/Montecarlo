from particles import Particle, Material
import numpy as np
from matplotlib import pyplot as plt

melec = 511

def bethe_bloch(p, medium):

    r0 = 2.82e-13 #[cm]
    Na = 6.022e23
    N = medium.density * Na / medium.A
    A = np.log(p.beta()*p.gamma()*np.sqrt(p.gamma()-1)*p.mass/medium.I)
    #B = (1 / (2*p.gamma())) * ((((p.gamma()-1)**2)/8) +1 - (2*p.gamma()**2 + 2*p.gamma()-1)*np.log(2))
    B = (1 - np.log(2)*(2*p.gamma()**2 + 2*p.gamma() - 1) + (p.gamma()**2 -1)/8   ) / (2*p.gamma()**2)

    let = 4*np.pi*(r0**2)*(p.mass/(p.beta()**2)) * N * medium.Z * (A + B)

    return let

def bethe_bloch2(p, medium):


    #C = 0.1536 # MeV cm ^2 / g
    C = 153.6 # keV cm ^2 / g
    tau = p.eKin / melec

    def F(t):

        if p.charge < 0:
            #print('Sono un elettrone')
            return 1 - p.beta()**2 + ((t**2/8) - (2*t + 1)*np.log(2)) / (t + 1)**2
        else:
            #print('Sono un positrone')
            return 2*np.log(2) - ((p.beta()**2)/12) * (23 + 14/(t+2) + 10/((t+2)**2) + 4/((t+2)**3))


    let = (C * medium.Z/medium.A * medium.density * p.charge**2 / p.beta()**2) * (np.log(  ((tau+2)*tau**2) / (2*(medium.I/melec)**2)  ) + F(tau))

    return let


if __name__ == '__main__':

    water = Material('Water', 7.2, 18, 1, 75e-3)

    p = Particle('Positron', 511, 1)
    p.eKin=np.linspace(0000.1, 510e3, 100000)

    plt.figure(1)
    plt.title('Collision energy loss; e$^+$ e$^-$ in $H_2O$')
    #plt.plot(e.eKin, bethe_bloch(e, water), label='e$^-$')
    #plt.plot(p.eKin, bethe_bloch1(p, water), label='1')
    plt.plot(p.eKin, bethe_bloch2(p, water), label='2')
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(0.1, 100)
    #plt.ylim(1.5, 10)
    plt.legend(loc='best')
    plt.xlabel('p / M c')
    plt.ylabel('Energy loss [MeV cm^2 / g]')
    plt.show()
import numpy as np
from matplotlib import pyplot as plt


LIGHT_SPEED = 1.


class Particle:
    ''' Class describing a particle'''

    def __init__(self, name, mass, charge, eKin=0):
        ''' Arguments:
            - name of the particle
            - mass [ MeV/c^2 ]
            - charge [ e ]
            - kinetic energy [MeV]
         '''
        self.name = name
        self.mass = mass
        self.charge = charge
        self.eKin = eKin

    def __str__(self):
        message = ' \n'
        message += 'Particle {}: \n'
        message += '- Mass = {} Mev/c^2 \n'
        message += '- Charge = {} e \n'
        message += '- Kinetic energy = {} MeV \n'
        return message.format(self.name, self.mass, self.charge, self.eKin)

    def energy(self):
        return self.eKin + (self.mass * LIGHT_SPEED**2)

    def gamma(self):
        return self.energy() / (self.mass * LIGHT_SPEED ** 2)

    def beta(self):
        return np.sqrt(1 - (1 / (self.gamma()**2)))

    def momentum(self):
        return np.sqrt(self.energy()**2 - (self.mass * LIGHT_SPEED**2)**2) / LIGHT_SPEED


class BetaSource:

    def __init__(self, name, Emax, Z):
        ''' Arguments:
            - name of the beta source
            - atomic number
            - maximum (end point) energy of the beta particles [keV]
         '''
        self.name = name
        self.Z = Z
        self.Emax = Emax


    def beta_spectrum(self, p):

        def fermi_function(p, Z_daughter):
            alpha = 1 / 137
            Z_daughter = self.Z - 1
            n = -Z_daughter * alpha * p.energy() / (p.momentum())
            return 2 * np.pi * n / (1 - np.exp(-2 * np.pi * n))

        w = 1 + p.eKin / p.mass
        w0 = 1 + self.Emax / p.mass
        momentum = np.sqrt(w**2 - 1)
        spectra = fermi_function(p, self.Z-1) * (p.momentum()/p.mass) * w * ((w0 - w)**2)

        return spectra


class Material:
    ''' Class describing a material'''

    def __init__(self, name, Z, A, density, I):
        ''' Arguments:
            - name of the material
            - atomic number
            - mass number
            - density [g / cm^3]
            - ionization energy [MeV]
        '''
        self.name = name
        self.Z = Z
        self.A = A
        self.density = density
        self.I = I

    def __str__(self):
        message = ' \n'
        message += 'Material {}: \n'
        message += '- Atomic number = {}  \n'
        message += '- Atomic weight = {}  \n'
        message += '- Density = {} (g / cm^3 ) \n'
        message += '- Ionization energy = {} MeV \n'
        return message.format(self.name, self.Z, self.A, self.density, self.I)


if __name__ == '__main__':

    p = Particle('Positron', 511, +1)

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.title('Theoretical energy spectra')
    plt.xlabel('Energy [keV]')
    plt.ylabel('Counts [a.u]')
    plt.ylim(0, 1.2)
    plt.grid()

    F18 = BetaSource('Fluoro', 635, 9)
    p.eKin = np.linspace(0.1, F18.Emax, 1000)
    plt.plot(p.eKin, F18.beta_spectrum(p)/max(F18.beta_spectrum(p)), color='red', label='$^{18}$F')

    C11 = BetaSource('Carbon', 970, 6)
    p.eKin = np.linspace(0.1, C11.Emax, 1000)
    plt.plot(p.eKin, C11.beta_spectrum(p)/max(C11.beta_spectrum(p)), color='lime', label='$^{11}$C')

    N13 = BetaSource('Azoto', 1190, 7)
    p.eKin = np.linspace(0.1, N13.Emax, 1000)
    plt.plot(p.eKin, N13.beta_spectrum(p)/max(N13.beta_spectrum(p)), color='blue', label='$^{13}$N')

    O15 = BetaSource('Ossigeno', 1720, 8)
    p.eKin = np.linspace(0.1, O15.Emax, 1000)
    plt.plot(p.eKin, O15.beta_spectrum(p)/max(O15.beta_spectrum(p)), color='cyan', label='$^{15}$O')

    plt.legend(loc='best')


    plt.show()
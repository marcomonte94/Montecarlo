import numpy as np

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
        return self.energy() / self.mass * LIGHT_SPEED ** 2

    def beta(self):
            return np.sqrt(1 - (1 / self.gamma()**2))




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
    p = Particle('Proton', 938., 1, 100)
    print(p.energy())
    water = Material('Water', 10, 18, 1, 75e-6)
    print(p)
    print(water)

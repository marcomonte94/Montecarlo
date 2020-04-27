import numpy as np 

class Particle:

    def __init__(self, name, mass, charge):
        self.name = name
        self.mass = mass
        self.charge = charge

    def __str__(self):
        message = ' \n'
        message += 'Particle {}: \n'
        message += '- Mass = {} Mev/c^2 \n'
        message += '- Charge = {} e \n'
        return message.format(self.name, self.mass, self.charge)

    def gamma(self, b):
        return 1 / np.sqrt(1 - b**2)

        
if __name__ == '__main__':
    e = Particle('electron', 0.511, 1)
    alpha = Particle('Alpha', 4*938, 2)
    print(e)
    print(alpha)

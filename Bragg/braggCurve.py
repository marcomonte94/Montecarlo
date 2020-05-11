import numpy as np
from particles import Particle, Material
from LinearEnergyTransfer import bethe_bloch
from matplotlib import pyplot as plt

melec = 0.511

np.random.seed(80)


def samplingTheta(p, mat, dx):
    ''' Scattering angle sampling following Lynch and Dahl's empirical formula
    '''
    '''
     IT DOESN'T WORK!!! WTF???

    F = 95/100
    xa2 = 2.007e-5 * mat.Z**(2/3) *(1 + 3.34*(mat.Z*1/(137*beta)**2))/E**2
    a = (mat.Z*(mat.Z+1)/mat.A)

    xc2 = 0.157 * 1 * a*(mat.density*dx/(E**2) * beta**2)
    omega = xc2 / xa2
    v = 0.5 * omega / (1 - F)

    thetaVar = 0.5*(2 * xc2 / (1+F)) * (((1+v)/v) * np.log(1+v)-1  )

    theta = np.random.normal(0, np.sqrt(thetaVar))
    '''

    ''' Scattering angle sampling following Molière's theory
    '''
    x0 = 36 # Radiation lenght [cm]
    theta0 = (13.6 / (p.beta()*p.eKin)) * 1 * np.sqrt(dx / x0)*(1 + 0.038*np.log(dx/x0))

    theta = np.random.normal(0, theta0)

    return theta


def enStraggling(p, mat, let, dx):
    ''' Gaussian energy straggling (only Gaussian for now...)
    '''

    Tmax = (2*melec*(p.beta()*p.gamma()**2))/ (1 + 2*p.beta()*(melec/p.mass)+(melec/p.mass)**2)
    k = let*dx / Tmax
    #print(k)

    sigma = 0.1569*mat.density * (mat.Z/mat.A)*dx
    sigma = sigma*(1-0.5*p.beta()**2)/(1-p.beta()**2)
    sigma = np.sqrt(sigma)
    de = -1
    while de < 0:
        de = np.random.normal(let*dx, sigma)

    return de


if __name__ == '__main__':

    water = Material('Water', 10, 18, 1, 75e-6)
    step = 0.005 #cm
    npart = 70 # n° of simulated particles
    r = []

    for i in range (npart):

        p = Particle('Proton', 938., 1, 100)

        ''' Initial energy Gaussian straggling
            (Mary Giusy docet)
        '''
        p.eKin = np.random.normal(p.eKin, 0.5)

        ''' Initial position Gaussian straggling
        '''
        x, y = np.random.normal(0, 0.001), np.random.normal(0, 0.001)
        #x, y = np.random.uniform(0, 0.001), np.random.uniform(0, 0.001)
        xpath, ypath = [], []
        eLoss = []

        ''' Multiple scattering until kinetic energy is over
        '''

        while p.eKin > 0:

            let = bethe_bloch(p, water)

            theta = samplingTheta(p, water, step)
            de = enStraggling(p, water, let, step)

            p.eKin = p.eKin - de
            x = x + step* np.cos(theta)
            y = y + step* np.sin(theta)

            eLoss.append(let)
            xpath.append(x)
            ypath.append(y)

        r.append(x)

        plt.figure(1)
        plt.title('Bragg curve: 70 proton in $H_20$, 100 MeV')
        plt.plot(xpath, eLoss, '.', color='blue')
        plt.ylim(0, 1000)
        plt.xlabel('x [cm]')
        plt.ylabel('Energy loss [MeV cm^2 / g]')

        plt.figure(2)
        plt.title('Proton path in xy plane')
        plt.plot(xpath, ypath, '.')
        plt.xlabel('x [cm]')
        plt.ylabel('y [cm]')

    plt.show()














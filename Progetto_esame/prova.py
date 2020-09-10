import numpy as np
from particles import Particle, BetaSource, Material
from eLoss import bethe_bloch, bremsstrahlung
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

melec = 511

def samplingTheta(p, mat, dx):
    
    
    '''
    F = 95/100
    xa2 = 2.007e-5 * mat.Z**(2/3) *(1 + 3.34*(mat.Z*1/(137*p.beta()**2)))/p.momentum()**2
    a = (mat.Z*(mat.Z+1))

    #xc2 = 0.157 * 1 * a*(mat.density*dx/(p.momentum()**2) * p.beta()**2)
    xc2 = 0.157 * p.charge * a * dx / ((p.momentum()**2) * (p.beta()**2))
    omega = xc2 / xa2
    v = 0.5 * omega / (1 - F)

    thetaVar = 0.5*(2 * xc2 / (1+F)) * (((1+v)/v) * np.log(1+v)-1  )
    theta = np.random.normal(0, np.sqrt(thetaVar))
    '''

    x0 = 36 # Radiation lenght [cm]
    theta0 = (13600 /( p.beta()*p.momentum())) * 1 * np.sqrt(dx / x0)*(1 + 0.038*np.log(dx/x0))

    theta = np.random.normal(0, theta0)
    
    return theta


def moliere(p, mat, dx):

    F = 95/100
    h = 6.58e-19 #[keV s]
    Na = 6.022e23
    a0 = 5.3e-8 #[cm]
    alpha = 1 / 137
    N = N = mat.density * Na / mat.A
    a = (mat.Z*(mat.Z+1))

    x0 = (h / p.momentum()) / (0.885*a0*mat.Z**(-1/3))
    xa2 = (1.13 + 3.76*(alpha**2))*(x0**2)
    xc2 = (4*np.pi*N*(np.e**4)*a*dx) / ((p.momentum()**2) * (p.beta()**2))

    

    return theta




if __name__ == '__main__':

    water = Material('Water', 7.22, 18, 1, 75e-3)
    p = Particle('Positron', 511, 1)
    p.eKin = 20e3
    theta = []

    for i in range(10000):
        theta.append(samplingTheta(p, water, 0.005))

    plt.hist(np.degrees(theta))
    plt.show()

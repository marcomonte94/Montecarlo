import numpy as np
from particles import Particle, BetaSource, Material
from eLoss import bremsstrahlung, bethe_bloch
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from scipy.optimize import curve_fit


melec = 511

def gaus(x, a, mu, sigma):
    return a * norm.pdf(x, mu, sigma)


def samplingEnergy(beta_sourge):

    p = Particle('Positron', 511, 1)
    i = 0

    while i != 1:
        p.eKin = np.random.uniform(0, beta_sourge.Emax) 
        y = np.random.uniform(0, 0.8)

        if y < beta_sourge.beta_spectrum(p):
            i = 1           

    return p.eKin

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
    theta0 = (13600 / p.eKin) * 1 * np.sqrt(dx / x0)*(1 + 0.038*np.log(dx/x0))

    theta = np.random.normal(0, theta0/np.sqrt(2))
    
    return theta

def enStraggling(p, mat, let, dx):
    ''' Gaussian energy straggling (only Gaussian for now...)
    '''

    Tmax = (2*melec*(p.beta()*p.gamma()**2))/ (1 + 2*p.beta()*(melec/p.mass)+(melec/p.mass)**2)
    k = let*dx / Tmax
    #print(k)

    sigma = 156.9*mat.density * (mat.Z/mat.A)*dx
    sigma = sigma*(1-0.5*p.beta()**2)/(1-p.beta()**2)
    sigma = np.sqrt(sigma)
    de = -1
    while de < 0:
        de = np.random.normal(let*dx, sigma)

    return de


if __name__ == '__main__':

    F18 = BetaSource('Fluoro', 635, 9)
    water = Material('Water', 7.22, 18, 1, 75e-3)
    Estep = 3
    ds = 0.1
    
    positrons = np.full(50000, Particle('Positron', 511, 1))
    '''
    sampled_spectrum = np.zeros(10000)

    for i in range(len(sampled_spectrum)):
        print(i)
        sampled_spectrum[i] = samplingEnergy(F18)

    plt.hist(sampled_spectrum, bins=50)
    plt.show()

    '''
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    x_endpoint, y_endpoint, z_endpoint = [], [], []
    #plt.figure()
    conta = 0

    for i_positron in positrons:
        conta += 1
        print(conta)
        xpath, ypath, zpath = [0], [0], [0]
        
        x, y, z = 0, 0, 0
        #print(i_positron.mass)
        i_positron.eKin = samplingEnergy(F18)
        #i_positron.eKin = 242
        r = np.random.uniform(-1, 1)
        theta = np.arccos(r)
        phi = np.random.uniform(0, 2*np.pi)
        #print(i_positron.eKin)

        while i_positron.eKin > 0.01:
            let = bethe_bloch(i_positron, water) #+ bremsstrahlung(i_positron, water)
            #ds = Estep / let
            theta += samplingTheta(i_positron, water, ds)
            #print(theta)
            phi += np.random.uniform(0, 2*np.pi)
            #de = enStraggling(i_positron, water, let, ds)
            de = let*ds
            i_positron.eKin = i_positron.eKin - de
            print(i_positron.eKin)

            x += ds * np.cos(theta) * np.sin(phi)
            y += ds * np.sin(theta) * np.sin(phi)
            z += ds * np.cos(phi)

            xpath.append(x)
            ypath.append(y)
            zpath.append(z)

        #x_endpoint.append(xpath[-1])
        #y_endpoint.append(ypath[-1])
        #z_endpoint.append(zpath[-1])

        x_endpoint.append(x)
        y_endpoint.append(y)
        z_endpoint.append(z)


 
        #ax.plot(xpath, ypath, zpath)
        plt.plot(xpath, ypath, '.', color='black')
        plt.plot(xpath, ypath, color='blue')

    plt.figure()
    plt.plot(y_endpoint, z_endpoint, '.')

    plt.figure()
    '''
    ydata, edges, _ = plt.hist(x_endpoint, bins=100, density='True')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    '''
    p0 = [max(ydata), 0, 1e-2]
    popt1, pcov1 = curve_fit(gaus, xdata, ydata)
    plt.plot(xdata, gaus(xdata, *popt1))

    plt.show()

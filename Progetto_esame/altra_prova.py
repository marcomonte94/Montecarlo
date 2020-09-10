import numpy as np
from particles import Particle, BetaSource, Material
from eLoss import bethe_bloch, bremsstrahlung
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm
from scipy.optimize import curve_fit

np.random.seed(3488563484)
melec = 511


def ms3d(alpha, beta, gamma, theta, phi):

    mu = np.cos(theta)

    if np.abs(gamma) != 1:
        a = np.sqrt((1-mu**2) / (1-gamma**2))
        alpha_sc = mu*alpha + a*(alpha*gamma*np.sin(phi) + beta*np.cos(phi))
        beta_sc = mu*beta + a*(beta*gamma*np.sin(phi) - alpha*np.cos(phi))
        gamma_sc = mu*gamma - a*np.sin(phi)*(1-gamma**2)

    else:
        b = np.sqrt(1-mu**2)
        alpha_sc = gamma*beta*np.cos(phi)
        beta_sc = b*np.sin(phi)
        gamma_sc = mu*gamma

    return np.array([alpha_sc, beta_sc, gamma_sc])



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
    Estep = 5
    
    positrons = np.full(10000, Particle('Positron', 511, 1))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
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

        r = np.random.uniform(-1, 1)
        theta = np.arccos(r)
        phi = np.random.uniform(0, 2*np.pi)
        alpha, beta, gamma = np.sin(phi) * np.sin(theta), np.cos(phi) * np.sin(theta), np.cos(theta)
        n = np.array([alpha, beta, gamma ])
        #print(i_positron.eKin)

        while i_positron.eKin > 25e-6:
            let = bethe_bloch(i_positron, water) + bremsstrahlung(i_positron, water)
            r = Estep / let
            theta_sc = samplingTheta(i_positron, water, r)
            #print(theta)
            phi_sc = np.random.uniform(0, 2*np.pi)
            n = ms3d(*n, theta_sc, phi_sc)
            
            #print(i_positron.eKin)

            x += r * n[0]
            y += r * n[1]
            z += r * n[2]

            xpath.append(x)
            ypath.append(y)
            zpath.append(z)

            de = enStraggling(i_positron, water, Estep/r, r)
            i_positron.eKin = i_positron.eKin - de
        #x_endpoint.append(xpath[-1])
        #y_endpoint.append(ypath[-1])
        #z_endpoint.append(zpath[-1])

        x_endpoint.append(x)
        y_endpoint.append(y)
        z_endpoint.append(z)


 
        ax.plot(xpath, ypath, zpath)
        #plt.plot(xpath, ypath, '.', color='black')
        #plt.plot(xpath, ypath, color='blue')


    plt.figure()
    plt.plot(x_endpoint, y_endpoint, '.')

    plt.figure()
    
    ydata, edges, _ = plt.hist(z_endpoint, bins=100)
    xdata = 0.5 * (edges[1:] + edges[:-1])
    '''
    p0 = [max(ydata), 0, 1e-2]
    popt1, pcov1 = curve_fit(gaus, xdata, ydata)
    plt.plot(xdata, gaus(xdata, *popt1))
    '''
    plt.show()

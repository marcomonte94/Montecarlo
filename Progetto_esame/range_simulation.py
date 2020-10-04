import numpy as np
from particles import Particle, BetaSource, Material
from eLoss import bethe_bloch, bremsstrahlung
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as mcolors


#np.random.seed(3488563484)
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



def positron_track(positronSource, mat, Estep):

    positron = Particle('Positron', 511, 1)  
    xpath, ypath, zpath = [0], [0], [0]
    x, y, z = 0, 0, 0
    positron.eKin = samplingEnergy(positronSource)

    r = np.random.uniform(-1, 1)
    theta = np.arccos(r)
    phi = np.random.uniform(0, 2*np.pi)
    alpha, beta, gamma = np.sin(phi) * np.sin(theta), np.cos(phi) * np.sin(theta), np.cos(theta)
    n = np.array([alpha, beta, gamma ])

    while positron.eKin > 25e-6:
        
        
        let = bethe_bloch(positron, mat) #+ bremsstrahlung(positron, mat)
        r = Estep / let
        theta_sc = samplingTheta(positron, mat, r)
        phi_sc = np.random.uniform(0, 2*np.pi)
        
        x += (r * n[0])
        y += (r * n[1])
        z += (r * n[2])
        
        n = ms3d(*n, theta_sc, phi_sc)
        
        
        
        xpath.append(x)
        ypath.append(y)
        zpath.append(z)

        #de = enStraggling(positron, mat, Estep/r, r)
        positron.eKin = positron.eKin - Estep

    return xpath, ypath, zpath
    



if __name__ == '__main__':

    F18 = BetaSource('Fluoro', 635, 9)
    O15 = BetaSource('Ossigeno', 1720, 8)
    C11 = BetaSource('Carbonio', 970, 6)
    N13 = BetaSource('Azoto', 1190, 7)


    water = Material('Water', 7.22, 18, 1, 75e-3)
    Estep = 20
    n_positron = 50000
    x_endpoint, y_endpoint, z_endpoint = [], [], []

    fig = plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    ax = fig.add_subplot(111, projection='3d')

    for i in range(n_positron):
        print(f'Positron n° {i}')
        xpath, ypath, zpath = positron_track(F18, water, Estep)
        x_endpoint.append(xpath[-1]*10) 
        y_endpoint.append(ypath[-1]*10) 
        z_endpoint.append(zpath[-1]*10)  

        xpath, ypath, zpath = np.array(xpath), np.array(ypath), np.array(zpath)
        ax.plot(xpath*10, ypath*10, zpath*10)
        
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_zlabel('z [mm]')

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.hist2d(x_endpoint, y_endpoint, bins=[80, 80], norm=mcolors.PowerNorm(0.3), cmap='hot')
    plt.xlabel('x [mm]')
    plt.ylabel('y [mm]')
    plt.colorbar()
    
    '''
    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.ylabel('Counts')
    plt.hist(x_endpoint, bins=100, color='blue')
    plt.xlabel('x-endpoint [mm]')
    
    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    plt.hist(y_endpoint, bins=100, color='blue')
    plt.ylabel('Counts')
    plt.xlabel('y-endpoint [mm]')
    
    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)   
    ydata, edges, _ = plt.hist(z_endpoint, bins=100, color='blue')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    plt.ylabel('Counts')
    plt.xlabel('z-endpoint [mm]')
    '''
    plt.show()

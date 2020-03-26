import numpy as np 
import time
import random
import math
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit

''' Campionamento di una gaussiana attraverso due metodi: il primo utilizza il teorema del limite centrale (TLC),
    grazie al quale si può affermare che la somma di n variabili appartenenti alla loro distribuzione si distribuisce
    come una gaussiana nel limite di grandi n; il secondo utilizza il metodo di inversione della cdf.

'''

def model(x, a, mu, sigma):
    ''' Modello (gaussiana) da fittare
    '''
    from scipy.stats import norm
    #return (a/((np.sqrt(2*np.pi) * sigma))) * np.exp(-(x - mu)**2 / 2*sigma**2)
    return a * norm.pdf(x, mu, sigma)

def tlc_sampling(size):
    ''' Campionamento che utilizza il TCL
    '''
    t0 = time.time()
    i, z = 0, []

    while i < size:
        z.append(sum(np.random.uniform(size=10)))
        i += 1
    
    plt.figure()
    y_data, edges, _ = plt.hist(z, bins=200)
    x_data = 0.5 * (edges[1:] + edges[:-1])

    popt, pcov = curve_fit(model, x_data, y_data, p0=[1000, 5, 1])
    print(popt)
    print(np.sqrt(pcov.diagonal()))
    mask = y_data > 0
    chi2 = sum(((y_data[mask] - model(x_data[mask], *popt)) / np.sqrt(y_data[mask]))**2.)

    plt.title('$\mu={:.2f}, \ \sigma={:.2f}, \ \chi^2/ndof={:.2f}$'.format(popt[1], popt[2], chi2/(mask.sum()-len(popt))))
    plt.plot(x_data, model(x_data, *popt), color='r')
    print('Elapses time: {:.3f}'.format(time.time() - t0))  
    print('Deviazione std attesa: {:.2f}'.format(math.sqrt(10/12)))

def right_sampling(size):
    ''' Campionamento con il metodo di inversione.
    '''
    t0 = time.time()
    i, pdf = 0, []
    '''
    while i < size:
        r = math.sqrt(- 2 * np.pi * math.log(random.random()))
        phi = 2 * math.pi * random.random()
        pdf.append(r * math.cos(phi))
        i += 1
    '''
    pdf = np.sqrt(- 2 * np.pi * np.log(np.random.uniform(size=size))) * np.cos(2 * np.pi * np.random.uniform(size=size))

    plt.figure()
    y_data, edges, _ = plt.hist(pdf, bins=200)
    x_data = 0.5 * (edges[1:] + edges[:-1])

    popt, pcov = curve_fit(model, x_data, y_data, p0=[1000, 0, 1])
    print(popt)
    print(np.sqrt(pcov.diagonal()))
    mask = y_data > 0
    chi2 = sum(((y_data[mask] - model(x_data[mask], *popt)) / np.sqrt(y_data[mask]))**2.)

    plt.title('$\mu={:.2f}, \ \sigma={:.2f}, \ \chi^2/ndof={:.2f}$'.format(popt[1], popt[2], chi2/(mask.sum()-len(popt))))
    plt.plot(x_data, model(x_data, *popt), color='r')
    print('Elapses time: {:.3f}'.format(time.time() - t0))

if __name__ == '__main__':
    n = 100000
    tlc_sampling(n)
    right_sampling(n)
    plt.show()
import numpy as np 
import random
import time
import math
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit


def sampling(size):
    ''' Campionamento con il metodo del rigetto
    '''
    t0 = time.time()

    def f(x):
        ''' Funzione da campioanre.
        '''
        return (math.exp(1 / (x ** 2 + 1))-1)/4.2

    def g(x):
        ''' Funzione maggiorante
        '''
        return 1 /(np.pi * (1 + x ** 2))

    conta = 0
    chiamate = 0
    xdata, pdf = [], []

    while conta < size:
        
        x = math.tan((math.pi/2) * (2*random.random()-1))
        y = random.uniform(0, g(x))
        chiamate +=1
        if y < f(x):
            xdata.append(x)
            pdf.append(y) 
            conta += 1
    print(chiamate)

    def model(x, a, b):
        return a * (np.exp(1 / (x ** 2 + b))-1)

    y_data, edges, _ = plt.hist(xdata, bins=150, range=(-5, 5), density=False)
    x_data = 0.5 * (edges[1:] + edges[:-1])
    popt, pcov = curve_fit(model, x_data, y_data)
    print('{:.2f}, {:.2f}'.format(*popt))
    print('{:.2f}, {:.2f}'.format(*np.sqrt(pcov.diagonal())))
    mask = y_data > 0
    chi2 = sum(((y_data[mask] - model(x_data[mask], *popt)) / np.sqrt(y_data[mask]))**2.)
    print(f'Chi quadro: {chi2/(len(x_data[mask])-len(popt))}')

    plt.title('Sampling. Elapses time: {:.3f} s'.format(time.time() - t0))
    _x = np.linspace(-5, 5, 1000)
    plt.plot(_x, model(_x, *popt), color='r')
    plt.xlabel('x')
    plt.ylabel('pdf(x)')

if __name__ == '__main__':
    n = 100000
    sampling(n)
    plt.show()
import numpy as np 
import time
import random
import math
from matplotlib import pyplot as plt 

def f(x):
    ''' Modello da fittare: esponenziale.
    '''
    return np.exp(-x)

def dumb_sampling(size):
    ''' Campionamento con il metodo del rigetto
    '''
    t0 = time.time()
    conta = 0
    chiamate = 0
    xdata, pdf = [], []

    while conta < size:
        x = 5 * random.random()
        y = random.random()
        chiamate +=1
        if y < f(x):
            xdata.append(x)
            pdf.append(y) 
            conta += 1
    
    plt.figure('Dumb sampling. Elapses time: {:.3f}'.format(time.time() - t0))
    plt.hist(xdata, bins=50, density=True)
    print(chiamate)

def right_sampling(size):
    ''' Campionamento con il metodo di inversione.
    '''
    t0 = time.time()
    pdf = []
    conta = 0
    while conta < size:
        y = - np.log(random.random())
        if y < 5:
            pdf.append(y)
            conta +=1
    plt.figure('Right sampling. Elapses time: {:.3f}'.format(time.time() - t0))
    plt.hist(pdf, bins=50, density=True)

def vectorized_sampling(size):
    ''' Campionamento con il metodo di inversione, vettorizzato
    '''
    t0 = time.time()
    pdf = - np.log(np.random.uniform(size=size))
    plt.figure('Right sampling. Elapses time: {:.6f}'.format(time.time() - t0))
    plt.hist(pdf, bins=50, density=True)
    

if __name__ == '__main__':
    random.seed(7)
    np.random.seed(7)
    dumb_sampling(100000)
    right_sampling(100000)
    vectorized_sampling(100000)
    plt.show()


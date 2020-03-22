import numpy as np 
import time
import random
import math
from matplotlib import pyplot as plt 

''' La funzione time.time() fornisce il numero di secondi 
    trascorsi dal 1° gennaio 1970, come nel C.
'''

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
        if y < math.exp(-x):
            xdata.append(x)
            pdf.append(y) 
            conta += 1
    
    plt.figure()
    plt.title('Dumb sampling. Elapses time: {:.3f} s'.format(time.time() - t0))
    plt.hist(xdata, bins=50, density=True)
    plt.xlabel('x')
    plt.ylabel('pdf(x)')
    print(f'Numero di chiamate alla funzione random.random(): {chiamate}')

def right_sampling(size):
    ''' Campionamento con il metodo di inversione.
    '''
    t0 = time.time()
    pdf = []
    conta = 0
    while conta < size:
        y = - math.log(random.random())
        if y < 5:
            pdf.append(y)
            conta +=1
    plt.figure()
    plt.title('Right sampling. Elapses time: {:.3f} s'.format(time.time() - t0))
    plt.hist(pdf, bins=50, density=True)
    plt.xlabel('x')
    plt.ylabel('pdf(x)')

def vectorized_sampling(size):
    ''' Campionamento con il metodo di inversione, vettorizzato
    '''
    t0 = time.time()
    pdf = - np.log(np.random.uniform(size=size))
    plt.figure()
    plt.title('Right sampling. Elapses time: {:.3f}'.format(time.time() - t0))
    plt.hist(pdf, bins=80, density=True)
    plt.xlabel('x')
    plt.ylabel('pdf(x)')
    plt.xlim(-0.5, 6)
    

if __name__ == '__main__':
    n = 100000
    random.seed(7)
    dumb_sampling(n)
    random.seed(7)
    right_sampling(n)
    np.random.seed(7)
    vectorized_sampling(n)
    plt.show()


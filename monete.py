import numpy as np 
import time
import argparse
from matplotlib import pyplot as plt 

_description='Programma di simulazione lanci di moneta.'

t0 = time.time()

def simulazione_moneta(lanci, ripetizioni):

    ris = []

    for i in range (ripetizioni):

        a = np.random.uniform(size=lanci)
        teste = a > 0.5
        ris.append(len(np.where(teste)[0]))  

    bin = np.random.binomial(lanci, 0.5, size=ripetizioni)

    plt.figure()
    plt.hist(ris, alpha=0.7, color='r', label='Risultati')
    plt.hist(bin, alpha=0.5, color='blue', label='Binomiale')
    plt.legend()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument('lanci', type=int, help='Lanci da effettuare')
    parser.add_argument('ripetizioni', type=int, help='Ripetizioni del numero di lanci')
    parser.add_argument('-g', '--graphic', help='Grafico istogrammi della simulazione e della previsione', action = 'store_true')
    args = parser.parse_args()
    simulazione_moneta(args.lanci, args.ripetizioni)

    print('Elapses time: {:.3f}'.format(time.time() - t0))

    if args.graphic:
        plt.show()

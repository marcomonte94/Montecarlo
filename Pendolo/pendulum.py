import numpy as np 
import random
import time
import math
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit

t0 = time.time() # start time
lgen = np.linspace(0.1, 1., 10)     # Pendulum length  (m)
g = 9.81 
thetamin = 0/57.3    # Minimum Theta in radiant
thetamax = 5/57.3 # Maximum Theta in radiant

class Student:
    
    def __init__(self, nmeas, accuracy, distract, reactivity):
        self.nmeas = nmeas #number of measures for each lenght
        self.accuracy = accuracy # student's accuracy
        self.distract = distract # student's distraction
        self.reactivity = reactivity # student's reactivity

def f(x, a):
    '''Fit model
    '''
    return a * x 

def no_isochronism(t):
    theta = theta = thetamin +(thetamax-thetamin)*np.random.uniform()
    t *= (1 + (theta ** 2) / 16)
    

def bad_counting(t, prob):
        r = random.random()
        if r < prob:
            t *= 0.9

def inerzia(t, l, r):
    k = 2/5
    t *= np.sqrt(1 + k * (r**2)/l**2)

def experiment(stud, inertia=False):
    '''Experiment workflow:

    1) Measuring lenght (lmis) and period (T)
    '''
    dl = 0.001
    lmis = lgen + 0.5 * dl * np.random.uniform(-1, 1, size=len(lgen))
    T, dT = np.zeros(len(lmis)), np.zeros(len(lmis))

    for i in range(len(lmis)):
        '''The student makes nmeas measures of the period
           for each lenght
        '''
        Tmeas = 2 * np.pi * np.sqrt(lmis[i] / g) + stud.reactivity * np.random.normal(size=stud.nmeas)

        if stud.accuracy==True:
            '''The student is accurated -> he considers the 
                contribute of non perfect isochronism.
            '''
            no_isochronism(Tmeas)

        if stud.distract:
            print(Tmeas)
            '''The student is distracted -> he makes
                wrong oscillation counting 
            '''
            bad_counting(Tmeas, 1)
            print(Tmeas)
            print('\n')

        if inertia:
            inerzia(Tmeas, lmis, 0.03)

        T[i] += Tmeas.mean()
        dT[i] += Tmeas.std()/np.sqrt(stud.nmeas)
        


    '''2) Recostruction data: fitting and chi-square calculation.
    '''
    T, dT = np.array(T), np.array(dT)
    xdata = lmis
    ydata = T**2
    d_xdata = dl 
    d_ydata = 2 * T * dT

    if stud.accuracy == False:
        '''The student isn't accurated -> he doesn't consider
           the uncertainty on lenghts
        '''
        xdata -= 0.01
        d_xdata = 0

    popt, pcov = curve_fit(f, xdata, ydata, sigma=d_ydata)
    grec = (4*np.pi**2) / popt
    d_grec = (grec/popt) * np.sqrt(pcov.diagonal())

    chi2 = sum(((ydata - f(xdata, *popt)) / np.sqrt(d_ydata**2 + (popt * d_xdata)**2))**2.)

    return grec, d_grec, chi2

def make_histo(occ, bins, results, description, color):
    '''Histogram building.
    '''
    plt.hist(occ, bins=bins, range=(0,100), alpha=0.5, color=color, label=description)
    plt.legend(loc='best')
    print('### Results for {}:'.format(description))
    print(" - Mean recostructed g: {:.2f} +- {:.2f}".format(results.mean(), results.std()))
    print(" - Mean chi square: {:.2f}".format(occ.mean()))
    print(" - STD chi square: {:.2f}".format(occ.std()))
    print('\n')


if __name__ == '__main__':
    t = 7
    no_isochronism(t)
    print(t)
    print('Eh, voleeevi!')
    

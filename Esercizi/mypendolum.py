import numpy as np 
import random
import time
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit

t0 = time.time() # start time
lgen = np.linspace(0.1, 1., 10)     # Pendulum length  (m)
g = 9.81 
thetamin = 5/57.3    # Minimum Theta in radiant
thetamax = 10/57.3 # Maximum Theta in radiant

class Student:
    
    def __init__(self, nmeas, accuracy, distract, tres):
        self.nmeas = nmeas #number of measures for each lenght
        self.accuracy = accuracy # student's accuracy
        self.distract = distract # student's distraction
        self.tres = tres # student's reactivity

def f(x, a):
    '''Fit model
    '''
    return a * x

def experiment(stud):
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
        Tmeas = 2 * np.pi * np.sqrt(lmis[i] / g) + stud.tres * np.random.normal(size=stud.nmeas)
        if stud.accuracy==True:
            '''The student is accurated -> he considers the 
                contribute of non perfect isochronism.
            '''
            theta = thetamin +(thetamax-thetamin)*np.random.uniform()
            Tmeas *= (1 + (theta**2)/16)
        if stud.distract:
            '''The student is distracted -> he makes
                wrong oscillation counting 
            '''
            Tmeas *= 0.9
        T[i] += Tmeas.mean()
        dT[i] += Tmeas.std()/np.sqrt(stud.nmeas)

    '''2) Recostruction data: fitting and chi-square calculation.
    '''
    T, dT = np.array(T), np.array(dT)
    xdata = lmis
    ydata = T**2
    d_ydata = 2 * T * dT
    popt, pcov = curve_fit(f, xdata, ydata, sigma=d_ydata)
    grec = (4*np.pi**2) / popt
    d_grec = (grec/popt) * np.sqrt(pcov.diagonal())

    if stud.accuracy == False:
        '''The student isn't accurated -> he doesn't consider
            the uncertainty on lenghts
        '''
        dl = 0
    chi2 = sum(((ydata - f(xdata, *popt)) / np.sqrt(d_ydata**2 + (popt * dl)**2))**2.)

    return grec, d_grec, chi2

def make_histo(occ, bins, results, description, color):
    '''Histogram building.
    '''
    plt.hist(occ, bins=bins, alpha=0.5, color=color, label=description)
    plt.legend(loc='best')
    print('### Results for {}:'.format(description))
    print(" - Mean recostructed g: {:.2f} +- {:.2f}".format(results.mean(), results.std()))
    print(" - Mean chi square: {:.2f}".format(occ.mean()))
    print(" - STD chi square: {:.2f}".format(occ.std()))
    print('\n')


if __name__ == '__main__':
    np.random.seed(3141456242)   # random seed
    nstud = 5000
    grec, chi2 = np.zeros(nstud), np.zeros(nstud)

    plt.figure()
    '''5000 ideal students simulation.
    '''
    for i in range(nstud):
        marco = Student(10, True, False, 0.01)
        grec[i] += experiment(marco)[0]
        chi2[i] += experiment(marco)[2]
    make_histo(chi2, 50, grec, 'Ideal student', 'blue')
    
    grec, chi2 = np.zeros(nstud), np.zeros(nstud)
    '''5000 lazy students simulation.
    '''
    for i in range(nstud):
        marco = Student(3, True, False, 0.01)
        grec[i] += experiment(marco)[0]
        chi2[i] += experiment(marco)[2]
    make_histo(chi2, 200, grec, 'Lazy student', 'red')

    grec, chi2 = np.zeros(nstud), np.zeros(nstud)
    '''5000 distracted students simulation.
    '''
    for i in range(nstud):
        marco = Student(10, True, True, 0.01)
        grec[i] += experiment(marco)[0]
        chi2[i] += experiment(marco)[2]
    make_histo(chi2, 50, grec, 'Distracted student', 'green')

    plt.xlim(0, 80)
    plt.xlabel('$\chi^2$')
    plt.ylabel('Occorrences')
    plt.legend(loc='best')
    print('Elapses time: {:.3f}'.format(time.time()-t0))
    plt.show()
    
        
    

    
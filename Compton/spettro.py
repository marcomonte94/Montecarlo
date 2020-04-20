import comptonFunction
import numpy as np
from setup import Detector
from matplotlib import pyplot as plt 

''' Grandezze (in metri) relative alle componenti
    dell'apparato sperimentale
'''
rcoll = 0.001 # raggio del collimatore
lcoll = 0.02 # lunghezza del collimatore
dcs = 0.5 # distanza collimatore - scintillatore
hscint = 1. #altezza scintillatore


myScintillator = Detector()
theta_ok = np.arctan(rcoll / lcoll)
nevents = 300000
i = 0
energy = []
e_cs = 662

while i < nevents:

    #Nascita di un fotone
    theta_start = np.random.uniform(0, np.pi / 2)
    pin = np.array([e_cs, e_cs*np.cos(theta_start), e_cs*np.sin(theta_start), 0])

    if np.abs(theta_start) < theta_ok:

        if (lcoll + dcs)*np.tan(theta_start) < hscint/2:

            p1 = pin
            pout = myScintillator.detect_photon(p1)


            trigger = myScintillator.check_momentum(pout) and myScintillator.check_energy(p1, pout)
            angle = np.arcsin(pout[2]/pout[0])
            angle = np.degrees(angle)

            if trigger:
                if (np.abs(angle-0)<20):
                    energy.append(myScintillator.acquireEnergy(pout))
    else:
        e1, psi = comptonFunction.SamplingKleinNishina(pin[0])
        theta_out = psi - np.pi - theta_start
        theta_max = np.arctan(2*rcoll / (lcoll - rcoll/np.tan(theta_start)))

        if theta_out >= 0 and theta_out < theta_max:
            p1 = comptonFunction.scatteringMatrix(pin, e1, psi)

            if np.tan(theta_out) < (hscint/2 + rcoll) / (lcoll + dcs - rcoll/np.tan(theta_start)):

                pout = myScintillator.detect_photon(p1)

                trigger = myScintillator.check_momentum(pout) and myScintillator.check_energy(p1, pout)
                angle = np.arcsin(pout[2]/pout[0])
                angle = np.degrees(angle)
                if trigger:
                    if (np.abs(angle-0)<40):
                        energy.append(myScintillator.acquireEnergy(pout))

    i += 1



plt.figure()
ydata, edges, _ = plt.hist(energy, bins=50, label='Energy distribution')
xdata = 0.5 * (edges[1:] + edges[:-1])
plt.errorbar(xdata, ydata, np.sqrt(ydata), fmt='.', color='black', label ='Data (Poissonian Error')
plt.legend(loc='best')
plt.xlabel('Energy [keV]')
plt.ylabel('Counts')
plt.show()











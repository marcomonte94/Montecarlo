import comptonFunction
import numpy as np
from setup import Detector
from matplotlib import pyplot as plt 

''' Grandezze (in metri) relative alle componenti
    dell'apparato sperimentale
'''
rcoll = 0.002 # raggio del collimatore
lcoll = 0.02 # lunghezza del collimatore
dcs = 1. # distanza collimatore - scintillatore
hscint = 0.05 #altezza scintillatore

np.random.seed(3488563484)
myScintillator = Detector()
theta_ok = np.arctan(rcoll / lcoll) # angolo massimo per non fare scattering nel collimatore
nevents = 300000
i = 0
angle = []
e_cs = 662

while i < nevents:

    #Nascita di un fotone
    theta_start = np.random.uniform(0, np.pi / 2)
    pin = np.array([e_cs, e_cs*np.cos(theta_start), e_cs*np.sin(theta_start), 0])

    if np.abs(theta_start) < theta_ok: # no scattering nel collimatore

        if (lcoll + dcs)*np.tan(theta_start) < hscint/2: # entra nel trigger?

            p1 = pin
            pout = myScintillator.detect_photon(p1)


            trigger = myScintillator.check_momentum(pout) and myScintillator.check_energy(p1, pout)

            if trigger:
                angle.append(myScintillator.acquireAngle(pout))

    else: # fa scattering nel collimatore
        e1, psi = comptonFunction.SamplingKleinNishina(pin[0])
        theta_out = psi - np.pi - theta_start # angolo di uscita rispetto all'asse x
        theta_max = np.arctan(2*rcoll / (lcoll - rcoll/np.tan(theta_start))) #angolo max per uscire

        if theta_out >= 0 and theta_out < theta_max:
            p1 = comptonFunction.scatteringMatrix(pin, eout, psi)

            if np.tan(theta_out) < (hscint/2 + rcoll) / (lcoll + dcs - rcoll/np.tan(theta_start)):
                # Entra nel trigger?
                pout = myScintillator.detect_photon(p1)

                trigger = myScintillator.check_momentum(pout) and myScintillator.check_energy(p1, pout)
                if trigger:
                    angle.append(myScintillator.acquireAngle(pout))


    i += 1

plt.figure()
print(len(angle))
ydata, edges, _ = plt.hist(angle, bins=50, label='Angular distribution')
xdata = 0.5 * (edges[1:] + edges[:-1])
plt.errorbar(xdata, ydata, np.sqrt(ydata), fmt='.', color='black', label ='Data (Poissonian Error')
plt.legend(loc='best')
plt.xlabel('Angle [°]')
plt.ylabel('Counts')
plt.show()











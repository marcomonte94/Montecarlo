import comptonFunction
import numpy as np
import matplotlib.pyplot as plt
import random

melec = 511

class Detector:
    '''Classe che descrive un rivelatore.
        Argomenti:
        - Soglia per trigger (theshold)
        - Risoluzione energetica
        - Risoluzione angolare
    '''
    def __init__(self, treshold=0, en_resolution=0.8, dtheta=5):
        ''' Costruttore
        '''
        self.treshold = treshold
        self.en_resolution = en_resolution
        self.dtheta = dtheta


    def check_momentum(self, pout):
        ''' Metodo che verifica che il fotone rivelato
            sia diretto in avanti (asse x)
        '''
        if pout[1] > 0:
            return True
        else: return False

    def check_energy(self, pin, pout):
        ''' Metodo che verifica che l'energia depositata
            dall'elettrone scatterato superi la soglia
            per poter dare il trigger
        '''
        pelec = np.array([melec,0,0,0])
        pelec = pin - pout + pelec
        theta = pelec[1]/np.sqrt(pelec[0]**2-melec**2)
        theta = np.arccos(theta)             #electron outgoing angle
        light = np.abs((pelec[0]-melec)/np.cos(theta))   # light output proportional to dE/dx
        r = np.random.normal()
        light = light + self.en_resolution * r * np.sqrt(light)
        if light > self.treshold:
            return True
        else: return False

    def detect_photon(self, p):
        ''' Metodo che misura il quadrimpulso del fotone scatterato
        '''
        eout, psi = comptonFunction.SamplingKleinNishina(p[0])
        pout = comptonFunction.scatteringMatrix(p, eout, psi)
        r = np.random.uniform()
        if r < 0.5:
            phi = 0
        else:
            phi = np.pi
        pout = comptonFunction.xRotation(pout, phi)
        return pout

    def acquireAngle(self, p):
        '''Metodo che calcola l'angolo (theta) che viene
            misurato dal secondo rivelatore.
        '''
        theta = p[2]/p[0]
        theta = np.arcsin(theta)
        theta = np.degrees(theta) #conversion in degree
        theta_meas = theta + np.random.uniform()*self.dtheta
        return theta_meas

    def acquireEnergy(self, p):
        ''' Metodo che calcola l'energia depositata dal fotone
            nel rivelatore, che serivrà per costruirne lo spettro
        '''
        r = np.random.uniform()
        ein = p[0]
        if (r>0.8): emeas = ein
        if (r<0.8):
            emax = 2*ein**2/(2*ein +melec)    #Compton shoulder 0<E<Emax
            r = np.random.uniform()
            emeas = r*emax
        r = np.random.normal()
        emeas = emeas+ r*self.en_resolution*np.sqrt(emeas)
        return emeas


if __name__ == '__main__':
    print('Definizione del rivelatore e delle sue caratteristiche')
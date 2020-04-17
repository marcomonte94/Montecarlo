import numpy as np

melec = 511

def SamplingKleinNishina(e0):
    '''Campionamento della sezione d'urto di K.N.
    '''
    epsilon = 2*e0/melec
    bingo = -1
    while bingo <0:
        r1,r2 = np.random.uniform(), np.random.uniform()
        x = 1+r1*epsilon              #x -box
        y = ((epsilon+2-2*x)/(epsilon*x))**2+1/x-1/x**2+1/x**3
        y = 0.5*y
        if r2<y:           #  Note that 0<y<1
            bingo=1
    a = 1-2*(x-1)/epsilon
    psi =np.arccos(a)
    return e0/x,psi

def scatteringMatrix(pin, eout, angle):
    '''Calcolo della matrice di scattering
    '''
    c, s = np.cos(angle), np.sin(angle)
    s_matrix = (eout/pin[0]) * np.array([[1,0,0,0],[0,c,-s,0],[0,s,c,0],[0,0,0,1]])
    return np.dot(s_matrix, pin)

def xRotation(pin, angle):
    '''Rotazione attorno all'asse x: qui decido se il fotone va nella parte superiore o inferiore
        del secondo rivelatore -> trasformo il range dell'angolo di scattering
        da (0°, 180°) a (-90°, 90°).
    '''
    c, s = np.cos(angle), np.sin(angle)
    r_matrix = np.array([[1,0,0,0],[0,1,0,0],[0,0,c,-s],[0,0,s,c]])
    return np.dot(r_matrix, pin)


def photonCompton(pin):
    '''Funzione che seleziona energia e angolo di scattering dalla distribuzione
        di K.N., e determina le componenti del quadrimpulso di uscita
        tramite la matrice di scattering.
    '''
    e_out, psi = SamplingKleinNishina(pin[0]) #energia e angolo di scattering
    p_out = scatteringMatrix(pin, e_out, psi) #4-impulso dopo scattering
    return p_out



if __name__ == '__main__':
    print('Funzioni da utilizzare per la simulazione dello scattering Compton.')

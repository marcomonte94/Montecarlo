import numpy as np
import math
import pandas as pd
import pylab as pl
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy.optimize import curve_fit
from datetime import datetime,timedelta
import matplotlib.pyplot as plt


#   1)  DEFINITION
#     Energy in Mev
melec = 511   #electron mass (KeV)
ecs   = 662   # Photon Energy (KeV) Source: Cesium 137
dtheta = 5    # angle resolution
de     = 0.5  # Energy/light resolution dE/E = 80%/sqrt(E)

#################################
#     Klein-Nishina Sampling
#     input : eo incident energy
#     output: eout output energy
#             psi =  scattering angle
######################################
def SamplingKleinNishina(e0):
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


#   2) Generation
##########################################

melec = 511             #electron mass (KeV)
ecs   = 662             # Photon Energy (KeV) Source: Cesium 137
dtheta = 5              # angle resolution
de     = 0.8            # Energy/light resolution dE/E = 80%/sqrt(E)
nsca = 16000              #total number of events
xout = np.zeros(nsca)
angle = np.zeros(nsca)
energy = np.zeros(nsca)

isca=0
while isca < nsca:
    trigger = 0           # trigger = 0 NOT/  YES >0
#    photon birth
    pin =([ecs,ecs,0,0])            #input photon 4-momentum
# compton Effect at Scintillator
    ein = pin[0]
    eout,psi= SamplingKleinNishina(ein) #Outgoing Photon Energy
#compute Scattering Matrix
    c,s = np.cos(psi),np.sin(psi)
    scatter = np.array([[1,0,0,0],[0,c,-s,0],[0,s,c,0],[0,0,0,1]])
    x = ein/eout
    scatter = scatter/x
# randomizing azimuthal angle
#################################################
# only left-right so far. This part IS to be modified for 3d-simulation
    r = np.random.uniform()
    if (r<0.5): phi=0
    if (r>0.5): phi = np.pi
#####################################################

# rotation around x -axis
    c,s =np.cos(phi),np.sin(phi)
    rot = np.array([[1,0,0,0],[0,1,0,0],[0,0,c,-s],[0,0,s,c]])
    pout = scatter.dot(pin)
    pout = rot.dot(pout)       #output photon 4-momentum
    if(pout[1]>0):trigger = trigger+1    #only forward photons accepted
#    tmp = pout[0]**2-pout[1]**2 -pout[2]**2-pout[3]**2
#    print('photon mass.....',tmp)

#   electron birth
    pelec = [melec,0,0,0]   #Initial Electron at rest
    pelec = pin - pout + pelec      #output electron momentum
#  Does the electron trigger ?  Simple model based on Light output
    theta = pelec[1]/np.sqrt(pelec[0]**2-melec**2)
    theta = np.arccos(theta)             #electron outgoing angle
    light = (pelec[0]-melec)/np.cos(theta)   # light output proportional to dE/dx
    r = np.random.normal()
    light = light + de*r*np.sqrt(light)     # dL/L = 80%/sqrt(L)  Resolution Model
    threshold = 0.                     # Set threshold on electron signal in Scintillator
    if (light>threshold): trigger = trigger+1
#   Photon angle Measured by NaI cristal  resolution +-5° (gaussian)
    theta = pout[2]/pout[0]
    theta = np.arcsin(theta)
    theta = theta*57.3 #conversion in degree
    r = np.random.uniform()
    theta_meas = theta + r*dtheta
    if (np.abs(theta-0)>30): trigger=0
#photon energy measured in NaI   resolution dE/E = +- 80%/sqrt(E)
#Here try to simulate Compton Shoulder effect
    r = np.random.uniform()
    ein = pout[0]
    if (r>0.8): emeas = ein
    if (r<0.8):
        emax = 2*ein**2/(2*ein +melec)    #Compton shoulder 0<E<Emax
        r = np.random.uniform()
        emeas = r*emax
    r = np.random.normal()
    emeas = emeas+ r*de*np.sqrt(emeas)

    if trigger==2:
         angle[isca] = theta_meas
         energy[isca] = emeas
         if (np.mod(isca,100)==0):print(100*isca/nsca,'%')
         isca = isca+1

#   3) Plotting  Angle Distribution
xout = energy
xmin = np.amin(xout)
xmax = np.amax(xout)
nbin =  50

# get the histos content (y and bins)
bins = np.linspace(xmin, xmax, nbin)
xfit = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
y, bins = np.histogram(xout,bins,range=(xmin,xmax), normed=None, weights=None, density=False)
y=y+1
dy = np.sqrt(y)   #Poissonian Error


#plot the histos
plt.errorbar(xfit,y,dy, fmt='.', color='black', label =r'Data (POissonian Error')
plt.hist(xout, bins,range=(xmin,xmax) ,density=False, weights =None,cumulative=False,color ='cyan',histtype='stepfilled',label = r' Simulation' )
print('mean = ', np.mean(xout))
print('standard deviation', np.std(xout))
#Make the plot nicer
plt.xlabel(r' x axis')
plt.ylabel(r'y-axis')
plt.legend(loc='best')
plt.title(' Compton Effect')
pl.xscale('linear')
pl.yscale('linear')


pl.show()
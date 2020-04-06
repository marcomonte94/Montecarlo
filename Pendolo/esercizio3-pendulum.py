import numpy as np
import scipy.stats as st
import matplotlib.pyplot as plt
from numpy import random
from scipy.optimize import curve_fit
import pylab as pl
import math
from sklearn.linear_model import LinearRegression
import sys

#  Simulation of LaB1 measurement of g
#   Using a simple Pendulum
#
#   10 Lenght
#   Student are divided in categories
#    Ideal:   Everything Fine
#    Lazy :   only Few measurements for each length   (nmeas = 3)
#    Inaccurate :  Wrong Error Propagation, No attention to small Oscillation, Bad Counting
#

# 1)  Initialization -  Initial value and other stuff
np.random.seed(3141456242)   # random seed
g = 9.81        #Real Value
nstu = 16000              #Number of students
chi2 = np.zeros(nstu)     # Foreach student a chi2 is returned
gvalue = np.zeros(nstu)   #for each student a  meausured g  is returned

nmeas= 10  #Measurements for each length
dt = 0.01   #Human resolution on reaction time (0.1 s Divided by 10!)(Seconds)
dl = 0.01   # Error on pendulum length
thetamin = 0/57.3    #Minimum Theta in radiant
thetamax = 5/57.3 #MAximum Theta in radiant
lgen = [.1,.2,.3,.4,.5,.6,.7,.8,.9,1.]     #Pendulum length  (m)

tgen = np.zeros((len(lgen),nmeas),dtype=float)
trec = np.zeros(len(lgen))
etrec= np.zeros(len(lgen))
lrec = np.zeros(len(lgen))
y =    np.zeros(len(lgen))
ey =   np.zeros(len(lgen))
x =    np.zeros(len(lgen))
ex =   np.zeros(len(lgen))

#    2) GENERATION PART
ist=0
while ist < nstu:
    for ilen in range(len(lgen)):
            r = np.random.uniform()
            r = 2*r-1
            lrec[ilen] = lgen[ilen]+0.5*dl*r       # smearing length
            for imeas in range(nmeas):
                a = 2*np.pi*np.sqrt(lgen[ilen]/g)
                theta = thetamin +(thetamax-thetamin)*np.random.uniform() #Initial angle smearing
                a = a*(1+theta*theta/16)    #Non perfect isochronism
                a = a+dt*np.random.normal() #gaussian smearing for resolution
                r = np.random.uniform()
                if r<0.66 : a = 0.9*a        #BAD counting simulation
                tgen[ilen][imeas] =a

# 3)    RECONSTRUCTION PART  (HERE SIMULATION OF WHAT STUDENT DO  )
#       HERE RECONSTRUCED VARIABLE ONLY CAN BE USED
    trec = np.mean(tgen, axis=1)
    etrec = np.std(tgen,axis=1)/np.sqrt(nmeas)     # error on mean
    thetarec = (thetamax-thetamin)/2.              #Correction for non isochronism
#    thetarec = 0.                                 #BAD No correction for non isochronism
    trec = trec/(1+thetarec*thetarec/16)           #REALLY GOOD BOY !! PROF. FEEL GOOD
# student collect data and call them x and y !
    x = lrec
#   x = lrec-0.01                                 # BAD not accounting for BOB size
#    dl=0                                          #BAD error propagation
    ex= np.ones(len(lgen), dtype=float)
    ex = dl*ex/np.sqrt(12)
    y,ey =  trec**2, 2*etrec*trec

    slope = y/x                              #  fitting y = slope*x     slope = (4*pi*pi/g)

    eslope = (ey/y)**2 + (ex/x)**2  #REALLY REALLY GOOD
    eslope = np.sqrt(eslope)
    eslope = eslope*slope
    a,b = sum(slope/(eslope*eslope)),sum(1/(eslope*eslope))          #weigthed  mean

    slope_rec = a/b                         #average slope
    eslope_rec = 1/np.sqrt(b)               #Error on slope

    a= (slope-slope_rec)/eslope
    chisq = sum(a*a)                        #Chi-square

    grec = 4*np.pi*np.pi
    grec= grec/slope_rec                   #Reconstructed g
    egrec = grec*eslope_rec/slope_rec      # Error on reconstructed g
    chi2[ist] = chisq                      #Store for plotting

    gvalue[ist] = grec
    ist = ist+1        #loop on student


#   4)    PLOTTING

xmin = 0
xmax = 100.
nbin = 100
# get the histos content (y and bins)
bins = np.linspace(xmin,xmax, nbin+1)
y, bins = np.histogram(chi2,bins,range=(xmin,xmax), normed=True, weights=None, density=None)
#plot the histos
plt.hist(chi2, bins,range=(xmin,xmax) ,density=True, weights =None,cumulative=False,color ='cyan',histtype='stepfilled',label = r' Simulated Chi-square' )
#superimpose the theoretical distribution (see chi-square distribution)
xfit = np.array([0.5 * (bins[i] + bins[i+1]) for i in range(len(bins)-1)])
ndof = len(lgen) -1   #number of degree of freedom

yfit =  xfit**(0.5*ndof -1)*np.exp(-0.5*xfit)
yfit = yfit/(2**(ndof/2))
yfit = yfit/math.gamma(ndof/2)
plt.plot(xfit, yfit, color='darkorange', linewidth=2.5, label=r' Expected chi-square')

# Make the plot nicer
plt.xlim(xmin,xmax)
ymax = np.amax(yfit)

pl.text(0.7*xmax,0.8*ymax, 'average %.3f'   % (np.mean(chi2)), fontsize=9)
pl.text(0.7*xmax,0.75*ymax,'st. dev %.3f' % (np.std(chi2)), fontsize=9)
pl.text(0.7*xmax,0.7*ymax,'g = %.3f +- %.3f'  % (np.mean(gvalue),np.std(gvalue)), fontsize=9)
plt.xlabel(r'Chi-square')
plt.ylabel(r'Arbitray Unit')
plt.legend(loc='best')
plt.title(' A g measurement with pendulum  ')

plt.show()
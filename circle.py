import numpy as np
import pandas as pd
import pylab as pl
import math
from numpy import exp,arange
from pylab import meshgrid,cm,imshow,contour,clabel,colorbar,axis,title,show
from scipy.optimize import curve_fit
from datetime import datetime,timedelta
from numpy import linalg as LA
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from scipy import special
from sklearn.linear_model import LinearRegression


#CIRCLE FITTING SIMULATION
#  Fitting Procedure
def circlefit(x,y):
   xb,yb = np.mean(x),np.mean(y)    # getting the baricenter
   u,v = x-xb,y-yb                   #passing to (u,v) coordinate system
   su,suu,suuu = sum(u),sum(u**2),sum(u**3)     #terms involving u
   sv,svv,svvv = sum(v),sum(v**2),sum(v**3)     #terms involving v
   suv= sum(u*v)                                 # second orders mixed terms
   suuv,suvv = sum(u*u*v),sum(u*v*v)             # third orders mixed terms

   uc = svv*(suuu+suvv)-suv*(svvv+suuv)
   uc = uc/(2*(suu*svv-suv**2))
   vc = suu*(svvv+suuv)-suv*(suuu+suvv)
   vc = vc/(2*(suu*svv-suv**2))

   rfit = uc**2+vc**2 + (suu+svv)/npoint
   xcfit,ycfit = xb+uc,yb+vc
   rfit = np.sqrt(rfit)                          #reconstructed Radius
   return  xcfit,ycfit,rfit


#  1)  Data simulation

ncircle = 100
diff = np.zeros(ncircle)
icircle = 0
while icircle < ncircle:
   xc,yc,r,ex,ey = -2,-2,10,0.5,0.5  # center coordinates x,y,radius, error on x, error on y
   npoint = 3   # Number of point
   phi = 2*np.pi*np.random.uniform(0,1,npoint)     # phi angle generation
   x,y = xc+r*np.cos(phi),yc+r*np.sin(phi)         # npoint coordinate generation
   dx,dy = np.random.uniform(-ex,ex,npoint),np.random.uniform(-ey,ey,npoint)
   x,y = x+dx,y+dy
   xcfit,ycfit,rfit = circlefit(x,y)
   diff[icircle] = r-rfit                           #Store true - reconstructed difference
   icircle = icircle+1

print('Media degli scarti = ', np.mean(diff))
rms= np.std(diff)
print('Deviazione Standard =' ,rms)



#) 1 plotting




#errorbar(x, y, yerr=None, xerr=None, fmt='', ecolor=None, elinewidth=None, capsize=None, barsabove=False, lolims=False, uplims=False, xlolims=False, xuplims=False, errorevery=1, capthick=None, *, data=None, **kwargs)
plt.xlim(-15,15)
plt.ylim(-15,15)
plt.errorbar(x,y,ex,ey, fmt='.', color='black', label =r'Generated point')
phi = t=np.linspace(0,2*np.pi,1000)
xfit,yfit = xc+rfit*np.cos(phi), yc+rfit*np.sin(phi)
plt.plot(xfit,yfit, color='darkorange', linewidth=2.5, linestyle='solid', label=r' Best Fitting circle')

# Make the plot nicer
xmax,ymax = np.amax(xfit),np.amax(yfit)
xmin,ymin = np.amin(xfit),np.amin(yfit)
#pl.text(0.7*xmax,0.7*ymax,'g = %.2f +- %.2f'  % (np.mean(gvalue),np.std(gvalue)), fontsize=9)
pl.text(3.8,9.5,'Radius =  %.3f +- %.3f ' %(rfit,rms), fontsize=9)
plt.xlabel(r' x axis')
plt.ylabel(r'y-axis')
plt.legend(loc='best')
plt.title(' Fitting a circle')
pl.xscale('linear')
pl.yscale('linear')

pl.show()


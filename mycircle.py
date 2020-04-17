import numpy as np
import pylab as pl
import math
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def circlefit(x,y, npoint):
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


def circle_simulation(xc, yc, r, npoint):

    ex,ey = 0.5, 0.5 # uncertainty
    phi = 2*np.pi*np.random.uniform(0, 1, npoint)     # phi angle generation
    x,y = xc + r * np.cos(phi), yc + r * np.sin(phi)  # npoint coordinate generation
    dx, dy = np.random.uniform(-ex ,ex, npoint), np.random.uniform(-ey, ey, npoint)
    x, y = x + dx, y + dy
    return circlefit(x, y, 10)  # Recostruction fit: xcfit, ycfit, rfit

if __name__ == '__main__':

    ncircle = 100
    xc, yc = -2,-2 # center's coordinates
    r = 10 # circle's radius
    diff = np.zeros(ncircle)

    for icircle in range(ncircle):  
        _x, _y, _r = circle_simulation(xc, yc, r, 10)
        diff[icircle] = r - _r 

    print(f'Raggio ricostruito {_r}')
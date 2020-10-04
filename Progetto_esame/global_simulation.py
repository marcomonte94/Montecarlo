import numpy as np
from particles import Particle, BetaSource, Material
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.optimize import curve_fit
from range_simulation import positron_track

F18 = BetaSource('Fluoro', 635, 9)
O15 = BetaSource('Ossigeno', 1720, 8)
C11 = BetaSource('Carbonio', 970, 6)
N13 = BetaSource('Azoto', 1190, 7)

water = Material('Water', 7.22, 18, 1, 75e-3)
Estep = 5
n_positron = 100000
x_endpoint, y_endpoint, z_endpoint = [], [], []
filedata = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_F_5.txt'

f = open(filedata,"w")
f.write('# x    y   z \n')


for i in range(n_positron):
    print(f'Positron n° {i}')
    xpath, ypath, zpath = positron_track(F18, water, Estep)
    x_endpoint.append(xpath[-1])
    y_endpoint.append(ypath[-1])
    z_endpoint.append(zpath[-1])

    f.write(f'{xpath[-1]} {ypath[-1]}  {zpath[-1]} \n')

f.close()
'''
plt.figure()
O_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_O_de.txt'
plt.figure()
x_endpoint, y_endpoint, z_endpoint = np.loadtxt(O_data, unpack=True)*10
x_endpoint = np.array(x_endpoint)
x_endpoint *= 10
ydata, edges, _ = plt.hist(x_endpoint, bins=400, density=False)
xdata = 0.5 * (edges[1:] + edges[:-1])
bin_size = (max(xdata)-min(xdata)) / 600

xfit, yfit = xdata[xdata>0.], ydata[xdata>0.]

ymax = max(ydata)

a = np.where(yfit<ymax/2)
x_2 = xfit[a[0][0]]
print(f'FWHM: {2*x_2} +/- {2*bin_size}')

b = np.where(yfit<ymax/10)
x_10 = xfit[b[0][0]]
print(f'FWTM: {2*x_10} +/- {2*bin_size}')
'''
plt.show()


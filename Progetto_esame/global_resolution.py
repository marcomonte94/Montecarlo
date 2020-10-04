import numpy as np
from matplotlib import pyplot as plt
from spatial_resolution import fit_and_compute
from scipy.stats import norm

w = 2 # [mm]
D = 200 # [mm]

plt.figure()

fwhm_w = w / 2
fwhm_D = 0.0022 * D # [mm]

F_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_o_de.txt'
x_endpoint, y_endpoint, z_endpoint = np.loadtxt(F_data, unpack=True)*10
ydata, edges, _ = plt.hist(x_endpoint, bins=500, density=False)
ynorm = ydata.sum()
plt.close()


xdata = 0.5 * (edges[1:] + edges[:-1])
fwhm_r, _= fit_and_compute(xdata, ydata, 6.5)
plt.close()

plt.figure(figsize=[7., 5.])
plt.rc('font', size=12)
#plt.xlim(-7.5, 7.5)

fwhm_tot = np.sqrt(fwhm_D**2 + fwhm_r**2 + fwhm_w**2)
x_fwtm = max(norm.pdf(xdata, 0, fwhm_tot / 2.35))/10
ydata = norm.pdf(xdata, 0, fwhm_tot / 2.35)
print(f'FWHM total: {fwhm_tot}')
x_10 = xdata[ydata < x_fwtm]
x_10 = x_10[x_10>0]
print(f'FWTM total: {2*x_10[0]}')


plt.hist(x_endpoint, bins=400, density=True, alpha=0.3, color='blue', label='Positorn Range')

plt.plot(xdata, norm.pdf(xdata, 0, fwhm_w / 2.35), color='orange', label='Detector Size')
plt.plot(xdata, norm.pdf(xdata, 0, fwhm_D / 2.35), color='green', label='Non-Collinearity')
plt.plot(xdata, norm.pdf(xdata, 0, fwhm_tot / 2.35), color='red', label='Total')
plt.xlabel('x [mm]')
plt.ylabel('Arbitraty Units')

results = 'FWHM = {:.2f} mm \n FWTM = {:.2f} mm'.format(fwhm_tot, 2*x_10[0])
plt.text(1.,1.5, results,bbox=dict(facecolor='none',edgecolor='black',boxstyle='square'))

plt.legend(loc='best')
plt.show()
from scipy.optimize import curve_fit
import numpy as np
from matplotlib import pyplot as plt


def fit_and_compute(xdata, ydata, fit_edge):


    def fitfunc(x, a, c, k1, k2):
        y = (c*np.exp(-k1*x) + (1-c)*np.exp(-k2*x))*a
        return y


    xfit, yfit = xdata[xdata>=0], ydata[xdata>=0]
    xfit, yfit = xfit[xfit<fit_edge], yfit[xfit<fit_edge]
    xfit, yfit = xfit[yfit>0], yfit[yfit>0]

    p0 = [ydata[0], 0.6, 15, 1.3]

    popt, pcov = curve_fit(fitfunc, xfit, yfit, p0, sigma=np.sqrt(yfit))
    chi2 = sum(((yfit - fitfunc(xfit, *popt)) / np.sqrt(yfit))**2.)
    print(f'Fit param: a = {popt[0]}, c = {popt[1]}, k1 = {popt[2]}, k2 = {popt[3]}')
    print(np.sqrt(pcov.diagonal()))
    ndof = len(xfit)-4
    print(f'Chi square norm: {chi2/ndof}')


    plt.figure()
    _x = np.linspace(0, max(xfit), 1500)
    plt.errorbar(xfit,yfit,np.sqrt(yfit), fmt='.', color='black', capsize=2, elinewidth=0.5, label =r'Data (Poissonian Error)')
    plt.plot(xfit, fitfunc(xfit, *popt), color='red', label='Fit Function')
    plt.yscale('log')
    plt.xlabel('x [mm]')
    plt.ylabel('Counts')
    plt.legend(loc='best')
    plt.grid()

    [a, c, k1, k2] = popt
    y_fwhm = a/2
    y_fwtm = a/10

    _x = np.arange(0, max(xfit), 1e-6)
    mask2 = fitfunc(_x, *popt) < y_fwhm
    mask10 = fitfunc(_x, *popt) < y_fwtm

    print('eccolo')
    print(2*_x[mask2][0])
    print(2*_x[mask10][0])

    return 2*_x[mask2][0], 2*_x[mask10][0]


if __name__ == '__main__':

    print('Fluoro \n')
    F_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_F.txt'
    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    x_endpoint, y_endpoint, z_endpoint = np.loadtxt(F_data, unpack=True)*10

    ydata, edges, _ = plt.hist(x_endpoint, bins=400, density=False, color='blue')
    plt.xlabel('x [mm]')
    plt.ylabel('Counts')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    x2_f, x10_f = fit_and_compute(xdata, ydata, 1.7)
    '''
    print('Carbonio \n')
    C_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_C_de.txt'

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    x_endpoint, y_endpoint, z_endpoint = np.loadtxt(C_data, unpack=True)*10
    ydata, edges, _ = plt.hist(x_endpoint, bins=600, density=False, color='blue')
    plt.xlabel('x [mm]')
    plt.ylabel('Counts')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    x2_c, x10_c = fit_and_compute(xdata, ydata, 2.7)

    print('Azoto \n')
    N_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_N_de.txt'

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    x_endpoint, y_endpoint, z_endpoint = np.loadtxt(N_data, unpack=True)*10
    ydata, edges, _ = plt.hist(x_endpoint, bins=500, density=False, color='blue')
    plt.xlabel('x [mm]')
    plt.ylabel('Counts')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    x2_n, x10_n = fit_and_compute(xdata, ydata, 3.7)

    print('Ossigeno \n')
    O_data = 'C:/Users/39348/Desktop/Montecarlo/Progetto_esame/hist_range_O_de.txt'

    plt.figure(figsize=[7., 5.])
    plt.rc('font', size=12)
    x_endpoint, y_endpoint, z_endpoint = np.loadtxt(O_data, unpack=True)*10
    ydata, edges, _ = plt.hist(x_endpoint, bins=500, density=False, color='blue')
    plt.xlabel('x [mm]')
    plt.ylabel('Counts')
    xdata = 0.5 * (edges[1:] + edges[:-1])
    x2_o, x10_o = fit_and_compute(xdata, ydata, 5.)




    '''


    plt.show()


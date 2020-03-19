import numpy as np 
import time
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit

a = -10 * np.log(np.random.uniform(size=1000000))

y_data, edges, _ = plt.hist(a, bins=50, density=False)
x_data = 0.5 * (edges[1:] + edges[:-1])

def model(x, a, b):
    return a*np.exp(-x/b)

popt, pcov = curve_fit(model, x_data, y_data, p0=[0.5, 0.5])
print(popt)
print(np.sqrt(pcov.diagonal()))

mask = y_data > 0
chi2 = sum(((y_data[mask] - model(x_data[mask], *popt)) / np.sqrt(y_data[mask]))**2.)

nu = mask.sum() - 3
sigma = np.sqrt(2 * nu)
print("{:.3f}, {:.3f}, {:.3f}".format(chi2, nu, sigma))

plt.plot(x_data, model(x_data, *popt), color='r')
plt.show()
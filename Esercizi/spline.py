import numpy as np
from matplotlib import pyplot as plt 
from scipy.interpolate import InterpolatedUnivariateSpline

x = np.linspace(0., 1., 10)
y = np.exp(3. * x) * np.sin(3. * x)

plt.plot(x, y, 'ro')

s1 = InterpolatedUnivariateSpline(x, y, k=3)
xs = np.linspace(0., 1., 1000)
plt.plot(xs, s1(xs))

plt.show()

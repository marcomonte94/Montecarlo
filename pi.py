import numpy as np 
from matplotlib import pyplot as plt 
import time 

t0 = time.time()
n = 100000

x = np.random.uniform(size=n)
y = np.random.uniform(size=n)

mask = x**2 + y**2 < 1

plt.plot(x[mask], y[mask], 'o')

print(f'Risultato: {4 * len(x[mask]) / n}')
print(f'Elapses time: {time.time() - t0}')

print(np.random.seed())

plt.show()


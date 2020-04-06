import pendulum
import numpy as np 
from matplotlib import pyplot as plt 


np.random.seed(3141456242)   # random seed
nstud = 5000
grec, chi2 = np.zeros(nstud), np.zeros(nstud)

plt.figure()
print('5000 ideal students simulation including momentum of inertia effect.')

for i in range(nstud):
    student = pendulum.Student(10, True, False, 0.01)
    grec[i] += pendulum.experiment(student)[0]
    chi2[i] += pendulum.experiment(student)[2]
pendulum.make_histo(chi2, 50, grec, 'Ideal student', 'blue')

plt.show()
import pendulum
import numpy as np 
from matplotlib import pyplot as plt 


np.random.seed(3141456242)   # random seed
nstud = 5000
grec, chi2 = np.zeros(nstud), np.zeros(nstud)

plt.figure()
'''5000 ideal students simulation.
'''
for i in range(nstud):
    student = pendulum.Student(10, True, False, 0.01)
    grec[i] += pendulum.experiment(student)[0]
    chi2[i] += pendulum.experiment(student)[2]
pendulum.make_histo(chi2, 50, grec, 'Ideal student', 'blue')

'''5000 inaccurate students simulation.
'''
grec, chi2 = np.zeros(nstud), np.zeros(nstud)
for i in range(nstud):
    student = pendulum.Student(10, False, False, 0.01)
    grec[i] += pendulum.experiment(student)[0]
    chi2[i] += pendulum.experiment(student)[2]
pendulum.make_histo(chi2, 50, grec, 'Inaccurate student', 'red')

'''5000 lazy students simulation.
'''
grec, chi2 = np.zeros(nstud), np.zeros(nstud)
for i in range(nstud):
    student = pendulum.Student(3, True, False, 0.01)
    grec[i] += pendulum.experiment(student)[0]
    chi2[i] += pendulum.experiment(student)[2]
pendulum.make_histo(chi2, 50, grec, 'Lazy student', 'green')

grec, chi2 = np.zeros(nstud), np.zeros(nstud)
'''5000 distract students simulation.
'''
for i in range(nstud):
    student = pendulum.Student(10, True, True, 0.01)
    grec[i] += pendulum.experiment(student)[0]
    chi2[i] += pendulum.experiment(student)[2]
pendulum.make_histo(chi2, 50, grec, 'Distract student', 'orange')

plt.show()
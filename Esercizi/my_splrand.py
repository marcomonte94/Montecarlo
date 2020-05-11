import numpy as np 
from matplotlib import pyplot as plt 
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.optimize import curve_fit



class ProbabilityDensityFunction(InterpolatedUnivariateSpline):

    '''Classe che descrive una funzione di densità di probabilità
    '''

    def __init__(self, x, y):
        # Richiamo il costruttore della classe madre
        InterpolatedUnivariateSpline.__init__(self, x, y)
        #Definisco la cdf (funzione cumulativa)
        y_cdf = np.array([self.integral(x[0], x_cdf) for x_cdf in x])
        self.cdf = InterpolatedUnivariateSpline(x, y_cdf)
        # Defiisco la ppf (percent-point function)
        x_ppf, i_ppf = np.unique(y_cdf, return_index=True) # se pdf non è monotona
        y_ppf = x[i_ppf]
        self.ppf = InterpolatedUnivariateSpline(x_ppf, y_ppf)


    def prob(self, x1, x2):
        '''Usa la cdf per ritornarmi la probabilità di trovare
        la variabile aleatoria tra x1 e x2
        '''
        return self.cdf(x2) - self.cdf(x1)

    def sampling(self, size=10000):
        '''Campiona la pdf con estrazioni random
        '''
        return self.ppf(np.random.uniform(size=size))

if __name__ == '__main__':
    x = np.linspace(0., 1., 101)
    y = 2. * x
    pdf = ProbabilityDensityFunction(x, y)

    plt.figure('PDF')
    plt.plot(x, pdf(x))
    plt.xlabel('x')
    plt.ylabel('pdf(x')

    plt.figure('cdf')
    plt.plot(x, pdf.cdf(x))
    plt.xlabel('x')
    plt.ylabel('cdf(x)')

    plt.figure('ppf')
    q = np.linspace(0., 1., 250)
    plt.plot(q, pdf.ppf(q))
    plt.xlabel('q')
    plt.ylabel('ppf(q)')

    plt.figure('Sampling')
    rnd = pdf.sampling(1000000)
    plt.hist(rnd, bins=200)

    plt.show()

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad, odeint
from scipy.signal import sawtooth

L = 1
def sawtooth_wave(x):
    return sawtooth(2*np.pi*x/(2*L))

def fourier_coeff(n_terms):
    def integrand_a0(x):
        return sawtooth_wave(x) / L
    def integrand_an(x, n):
        return sawtooth_wave(x) * np.cos(n*np.pi*x/L)/L
    def integrand_bn(x, n):
        return sawtooth_wave(x) * np.sin(n*np.pi*x/L)/L

    a0,_ = quad(integrand_a0, -L, L)

    an = np.zeros(n_terms)
    bn = np.zeros(n_terms)
    for n in range(1, n_terms + 1):
        an[n-1],_ = quad(integrand_an, -L, L, args=(n, ))
        bn[n-1],_ = quad(integrand_bn, -L, L, args=(n, ))
    return a0, an, bn

n_terms = 10

a0, an, bn = fourier_coeff(n_terms)

x = np.linspace(-L,L, 1000)

def f_approx(a0, an, bn, x, n_terms):
    result = a0/2
    for n in range(1, n_terms+1):
        result += an[n-1]*np.cos(n*np.pi*x/L) + bn[n-1]*np.sin(n*np.pi*x/L)
    return result
plt.plot(x, f_approx(a0, an, bn, x, n_terms))
plt.show()

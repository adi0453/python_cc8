import numpy as np
from scipy.integrate import odeint, quad
import matplotlib.pyplot as plt
from scipy.signal import square

#A
def CR(q, t):
    return (C*V0  - q)/ R*C

R, C, V0 = 1, 1, 10
q0 = 0
t = np.linspace(0, 5, 100)

sol = odeint(CR, q0, t)
#plt.plot(t, sol)
#plt.show()


#B
L = np.pi
x = np.linspace(-np.pi, np.pi, 1000)
def square_wave(x):
    return square(2* np.pi * x/(2*L))

def fourier_coeff(n_terms):
    def integrand_a0(x):
        return square_wave(x) / L

    def integrand_an(x, n):
        return square_wave(x)*np.cos(n*np.pi*x/L)/L
    
    def integrand_bn(x, n):
        return square_wave(x)*np.sin(n*np.pi*x/L)/L

    a0,_ = quad(integrand_a0, -L, L)

    an = np.zeros(n_terms)
    bn = np.zeros(n_terms)

    for n in range(1, n_terms + 1):
        an[n-1],_ = quad(integrand_an, -L, L, args = (n, ))
        bn[n-1],_ = quad(integrand_bn, -L, L, args = (n, ))
    return a0, an, bn

n_terms = 10
a0, an, bn = fourier_coeff(n_terms)

def fourier_series(a0, an, bn, x, n_terms):
    sol = a0/2
    for n in range(1, n_terms+1):
        sol += an[n-1]*np.cos(n*np.pi*x/L) + bn[n-1]*np.sin(n*np.pi*x/L)
    return sol

f = square_wave(x)
f_approx = fourier_series(a0, an, bn, x, n_terms)
plt.plot(x, f)
plt.plot(x, f_approx)
plt.show()

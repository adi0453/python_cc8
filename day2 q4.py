import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

#A
def Gaussian(x, sig):
    return (1/2*sig*np.sqrt(2*np.pi)*np.exp(-(x-x0)**2/(2*sig**2)))
def integrand(x, sig):
    return (x**3-4*x+2)*Gaussian(x, sig)

sig_values = [0.1, 0.05, 0.001]
x0 = 1
integrals = []
x = np.linspace(0, 2, 1000)

for value in sig_values:
    plt.plot(x, Gaussian(x, value))

#plt.show()
for value in sig_values:
    sol,_ = quad(integrand, 0, np.pi, args=(value, ))
    print(sol)


#B
def f(t, name):
    #print(name)
    return np.exp(-t**2/2)
    

sol,_ = quad(f, 0, np.inf, args=("adarsh",))
print(sol)

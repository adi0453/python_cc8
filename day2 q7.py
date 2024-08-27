import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def f(x):
    return np.log(x)/x**2

x = np.linspace(1, 100, 200)
integral_values = []

for n in range(200):
    integral,_ = quad(f, 1, np.inf)
    integral_values.append(integral)
    
plt.plot(x, f(x))
plt.show()

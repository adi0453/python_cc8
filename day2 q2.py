import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.integrate import quad

#A
def func(u, t):
    theta, z = u
    f1 = z
    f2 = -(g / L)*np.sin(theta)
    return [f1, f2]

L, g = 10, 980
u0 = [np.pi/6, 0]
t = np.linspace(0, 5, 1000)

sol = odeint(func, u0, t)
theta = sol[:, -1]
print(theta)
plt.xlabel("t")
plt.ylabel("theta")
#plt.plot(t, theta)
plt.title("Graph for swinging pendulum")
#plt.show()
    
#B
def f(x):
    return x**2*np.exp(-x**2)

N = np.linspace(0.1, 1, 100)
integral_array = []

for n in N: 
    integral,_= quad(f, -n, n)
    integral_array.append(integral)
plt.plot(N, integral_array)
plt.show()


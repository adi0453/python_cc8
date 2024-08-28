import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

def f(u, x):
    y, z = u
    f1 = z
    f2 = -y * 2*z**2/y
    return [f1, f2]

x = np.linspace(-1, 1, 50)
y0 = 1/(np.e + 1/np.e)
y0 = y1
z1, z2 = 0.1, 0.3 #guess values

#so we take initial guess and solve ODE, how can we implement this in eulers method and what is w1 equivalent in eulers method

u0 = [y0, z1]
sol = odeint(f, u0, x)
w1 = sol[:, 0][-1]

tol = 0.0001

for i in range(1000):
    u0 = [y0, z2]
    sol = odeint(f, u0, x)
    w2 = sol[:, 0][-1]

    if abs(w2-w1) < tol:
        break
    z1, z2 = z2, z2 + (z2 - z1) * (y1 - w2)/(w2 - w1)
    w1 = w2

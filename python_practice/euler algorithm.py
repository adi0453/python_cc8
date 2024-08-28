import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

'''
y = y + h * f(x, y, z)
'''

# solve ODE using Euler algorithm
def f1(x, y, z):
    return z

def f2(x, y, z):
    return -lam*z - k*y

x, y, z = 0, 0, 1.0
lam, k = 0.5, 2.0
h = 0.1

xlist, ylist = [x], [y]

while x <= 35:
    y += h * f1(x, y, z)
    z += h * f2(x, y, z)
    x += h
    xlist.append(x)
    ylist.append(y)
print(ylist)
plt.subplot(221)
plt.plot(xlist, ylist)



# solve ODE using odeint()
def f(u, x):
    y, z = u
    f1 = z
    f2 = -lam*z - k*y
    return [f1, f2]

u0 = [0, 1.0]
lam, k = 0.5, 2.0
x = np.linspace(0, 200, 100)

sol = odeint(f, u0, x)
print(sol)
plt.subplot(222)
plt.plot(x, sol[:,0])
#plt.show()


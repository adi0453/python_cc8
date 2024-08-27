import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jn

'''
#A
k = 0.05
nx = 100
nt = 300
x = np.linspace(0, np.pi, nx)
dx = x[1] - x[0]
dt = 0.01

u = np.zeros((nx, nt))

u[0, :] = 0
u[-1, :] = 0
u[:, 0] = np.sin(x)

for n in range(1, nt - 1):
    for i in range(1, nx - 1):
        u[i, n+1] = u[i, n] + k*(dt/dx)**2 * (u[i+1, n] - 2*u[i, n] + u[i -1, n])


plt.plot(x, u[:, 0], label='t=0')
plt.plot(x, u[:, int(nt/4)], label='t=0.75')
plt.plot(x, u[:, int(nt/2)], label='t=1.5')
plt.plot(x, u[:, -1], label='t=3.0')

plt.xlabel('x')
plt.ylabel('U(x,t)')
plt.legend()
plt.title('1D Heat Equation')
plt.show()
'''

#B
x = np.linspace(0, 10, 100)
n_values = [1, 2]

for n in n_values:
    Jn = jn(-n, x)
    Jn_negative = jn(-n, x)
    plt.plot(x, Jn)
    plt.plot(x, -Jn_negative)
plt.show()

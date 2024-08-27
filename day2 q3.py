import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

'''
#A
x = np.linspace(0, 1, 100)
y = np.linspace(0, 1, 100)

un = np.zeros((100, 100))


#boundary conditions
un[:, 0] = 0.0
un[:, -1] = y
un[0, :] = x
un[-1, :] = 0.0

for n in range(100):
    for i in range(1, 99):
        for j in range(1, 99):
            un[i, j] = 0.25*(un[i+1,j] + un[i-1, j] + un[i, j+1] + un[i, j-1])


def exact_sol(x, y):
    ua = x*y
    return ua
#plt.plot(x, un[50, :])
#plt.show()
'''
#B
def integrand(x):
    return x**2/np.sqrt((x-3)*(5-x))

sol,_ = quad(integrand, 3, 5)
x = np.linspace(3.001, 4.999, 400)
y = integrand(x)
print(sol)
plt.plot(x, y)
plt.show()

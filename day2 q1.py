from scipy.special import jn, j0
import numpy as np
import matplotlib.pyplot as plt

#A
x = np.linspace(1, 20,100)
J_0 = jn(0, x)
J_1 = jn(1, x)

def recurrence_relation(n, x, J_n, J_n_minus_1):
    return (2*n/x)*J_n - J_n_minus_1

J_2 = recurrence_relation(1, x, J_1, J_0)
J_3 = recurrence_relation(2, x, J_2, J_1)

#plt.plot(x, J_2)
#plt.plot(x, J_3)
#plt.show()


#B

x = np.linspace(0, 10, 101)
u = np.zeros(101)

u[100] = 100

for j in range(100):
    for i in range(100):
        u[i] += (u[i + 1] - 2*u[i] + u[i - 1])/4

plt.plot(x, u)
plt.show()

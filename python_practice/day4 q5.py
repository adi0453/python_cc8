import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

n1s = 1
n2s = 2
l = 0

def V(r):
    return (l*(l+1))/r**2 - 2 / r

r = np.linspace(0.1, 20, 1000)
dr = r[1] - r[0]
vr = V(r)

def eq_solver(E):
    u = np.zeros_like(r)
    u[0] = 1e-5
    u[1] = 1e-5

    for i in range(1, len(r) - 1):
        u[i+1] = 2*u[i] - u[i-1] + dr**2 * (vr[i] - E)*u[i]

    return u

def shooting(n):
    eigenvalues = []
    guess_energy = np.linspace(0.1, 1, 50)

    for E in guess_energy:
        u = eq_solver(E)
        if np.abs(u[-1]) < 1e-2:
            eigenvalues.append(E)
            if len(eigenvalues) >=  2.0:
                break
    return eigenvalues

def normalizer(u):
    norm = np.trapz(u**2 * r**2, r)
    return u / np.sqrt(norm)

eigenvalues_1s = shooting(n1s)[0]
eigenvalues_2s = shooting(n2s)[0]
print(eigenvalues_1s, eigenvalues_2s)

u1s = eq_solver(eigenvalues_1s)
u2s = eq_solver(eigenvalues_2s)

u1s_normalized = normalizer(u1s)
u2s_normalized = normalizer(u2s)

result = np.trapz(u1s_normalized * u2s_normalized * r**2, r)
print(result)

'''
for E in eigenvalues:
    u = eq_solver(E)
    u_normalized = normalizer(u)
    plt.plot(r, u_normalized)
plt.show()
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

m = 1
hbar = 1

def V(x):
    if x >= 0:
        return x
    else:
        return np.inf

x = np.linspace(0, 5, 1000)
vx = np.array([V(xi)for xi in x])
dx = x[1] - x[0]

def schrodinger_equation_solver(E):
    psi = np.zeros_like(x)
    psi[0] = 0
    psi[1] = 1e-5

    for i in range(1, len(x)-1):
        psi[i+1] = 2*psi[i] - psi[i-1] + dx**2 * (vx[i] - E)*psi[i]
    return psi

def shooting_method():
    eigenvalues = []
    guess_energies = np.linspace(0.5, 10, 100)

    for E in guess_energies:
        psi = schrodinger_equation_solver(E)
        if abs(psi[-1]) < 1e-2:
            eigenvalues.append(E)
            if len(eigenvalues) >= 2.0:
                break
    return eigenvalues

eigenvalues = shooting_method()

def normalizer(psi):
    norm,_ = quad(lambda xi: np.interp(xi, x, psi)**2, 0, 5)
    return psi / np.sqrt(norm)

for E in eigenvalues:
    psi = schrodinger_equation_solver(E)
    normalized = normalizer(psi)
    plt.plot(x, normalized)
plt.show()

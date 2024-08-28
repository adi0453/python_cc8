import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

m = 1
hbar = 1
omega = 1

def V(x):
    return 0.5 * m * omega**2 * x**2

x = np.linspace(-5, 5, 1000)
dx = x[1] - x[0]
vx = V(x)

def se_solver(E):
    psi = np.zeros_like(x)
    psi[0] = 1e-5
    psi[1] = 1e-5

    for i in range(1, len(x)-1):
        psi[i+1] = 2*psi[i] - psi[i-1] + dx**2 * 2*m/hbar**2 * (vx[i] - E) * psi[i]
    return psi

def shooting_function():
    eigenvalues = []
    guess_energies = np.linspace(0.5, 5, 100)

    for E in guess_energies:
        psi = se_solver(E)

        if np.abs(psi[-1]) < 1e-2:
            eigenvalues.append(E)
        if len(eigenvalues) >= 2.0:
            break
    return eigenvalues

def psi_normalizer(psi):
    N,_ = quad(lambda xi: np.interp(xi, x, psi)**2, -5, 5)
    return psi/np.sqrt(N)
eigenvalues = shooting_function()

for E in eigenvalues:
    psi = se_solver(E)
    psi_normalized = psi_normalizer(psi)
    plt.plot(x, psi_normalized)

plt.show()
    

import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

m = 1
hbar = 1
V0 = 40
a = 1

def V(x):
    if 0<= x <=a:
        return -V0
    else:
        return 0

x = np.linspace(-a, a, 1000)
dx = x[1] - x[0]
vx = np.array([V(xi) for xi in x])

def schrodinger_equation_solver(E):
    psi = np.zeros_like(x) #form an array like x
    psi[0] = 1e-5 #very small initial value
    psi[1] = 1e-5 #very small initial value

    #implementation of euler's algorithm
    for i in range(1, len(x)-1):
        psi[i+1] = 2*psi[i] - psi[i-1] + dx**2*(2*m/hbar**2)*(vx[i] - E)*psi[i]
    return psi


def shooting_algorithm():
    eigenvalues = []
    guess_energies = np.linspace(-45, -35, 100)

    for E in guess_energies:
        psi = schrodinger_equation_solver(E)
        if np.abs(psi[-1]) < 1e-2:
            eigenvalues.append(E)
            if len(eigenvalues) >= 2.0:
                break
    return sorted(eigenvalues)[:2]
            
eigenvalues = shooting_algorithm()

for E in eigenvalues:
    psi = schrodinger_equation_solver(E)
    norm_factor = np.sqrt(quad(lambda xi: np.interp(xi,x, psi,)**2, -a, a)[0])
    psi_normalized = psi / norm_factor
    plt.plot(x, psi_normalized, label=f'E = {E:.2f}')
plt.show()




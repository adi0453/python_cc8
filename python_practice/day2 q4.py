import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

D = 37244/10.5930
alpha = 2.380

def V(x):
    return D*(-1 + alpha**2 * x**2)

x_max = 0.5
x = np.linspace(-x_max, x_max, 1000)
vx = V(x)
dx = x[1] - x[0]

def schrodinger_equation_solver(E):
    psi = np.zeros_like(x)

    psi[0] = 1e-5
    psi[1] = 1e-5

    for i in range(1, len(x)-1):
        psi[i+1] = 2*psi[i] - psi[i-1] + dx**2 * (vx[i] - E)*psi[i]

    return psi

def shooting():
    eigenvalues = []
    guess_energies = np.linspace(0.1, 0.5, 100)

    for E in guess_energies:
        psi = schrodinger_equation_solver(E)

        if abs(psi[-1]) < 1:
            eigenvalues.append(E)
            if len(eigenvalues) >= 2.0:
                break
    return eigenvalues


def normalizer(psi):
    norm,_ = quad(lambda xi: np.interp(xi, x, psi)**2, -x_max, x_max)
    return psi / np.sqrt(norm)

eigenvalues = shooting()
print(eigenvalues)

for E in eigenvalues:
    psi = schrodinger_equation_solver(E)
    normalized_psi = normalizer(psi)
    plt.plot(x, normalized_psi)
    
plt.show()



    

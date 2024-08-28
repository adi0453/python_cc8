import numpy as np
import matplotlib.pyplot as plt

n = 1
l = 0
En = 1/n**2

def V(r):
    return (l*(l+1))/r**2 - 2/r

r_max = 20
r = np.linspace(0.1, r_max, 1000)
dr = r[1] - r[0]
vr = V(r)

def radial_equation_solver(E):
    u = np.zeros_like(r)
    u[0] = 1e-5
    u[1] = 1e-5

    for i in range(1, len(r)-1):
        u[i+1] = 2*u[i] - u[i-1] + dr**2 * (vr[i] - E) * u[i]

    return u

def shooting_method():
    eigenvalues = []
    guess_energies = np.linspace(0.1, 1.5, 50)

    for E in guess_energies:
        u = radial_equation_solver(E)

        if np.abs(u[-1]) < 1e-2:
            eigenvalues.append(E)
            if len(eigenvalues) >= 2.0:
                break
    return eigenvalues

eigenvalues = shooting_method()
print(eigenvalues)

def u_normalizer(u):
    integral = np.trapz(u**2 * r**2, r)
    return u / np.sqrt(integral)


for E in eigenvalues:
    u = radial_equation_solver(E)
    u_normalized = u_normalizer(u)
    plt.plot(r, u_normalized, label=f"E = {E:.2f}")
    plt.plot(r, u_normalized**2, label=f"probability density for E = {E:.2f}")
plt.legend()
plt.show()



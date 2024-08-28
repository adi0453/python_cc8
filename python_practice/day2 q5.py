import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Constants
l = 0  # Angular momentum quantum number for 2s orbital
n = 2  # Principal quantum number for 2s orbital

# Potential function
def V(r):
    return (l * (l + 1)) / r**2 - 2 / r

# Dimensionless radial Schrödinger equation solver using Euler's method
def schrodinger_equation_solver(E):
    u = np.zeros_like(r)
    u[0] = 1e-5  # Initial condition
    u[1] = 1e-5  # Small step to avoid singularity

    for i in range(1, len(r) - 1):
        u[i + 1] = 2 * u[i] - u[i - 1] + dx**2 * (V(r[i]) - E) * u[i]

    return u

# Shooting algorithm to find eigenvalues
def shooting_algorithm():
    eigenvalues = []
    guess_energies = np.linspace(0.1, 1.5, 50)

    for E in guess_energies:
        u = schrodinger_equation_solver(E)
        
        # Check if u(r) is close to zero at the boundary (r_max)
        if np.abs(u[-1]) < 1e-2:
            eigenvalues.append(E)
            if len(eigenvalues) >= 2:
                break

    return eigenvalues

# Normalize the radial wave function
def normalizer(u):
    # Integrate u(r)^2 * r^2 using the trapezoidal rule
    integral = np.trapz(u**2 * r**2, r)
    return u / np.sqrt(integral)

# Define the grid
r_max = 20
r = np.linspace(0.1, r_max, 1000)  # Avoid r = 0 to prevent division by zero
dx = r[1] - r[0]

# Find eigenvalues
eigenvalues = shooting_algorithm()
print("Eigenvalues:", eigenvalues)

# Plotting
for E in eigenvalues:
    u = schrodinger_equation_solver(E)
    u_normalized = normalizer(u)
    plt.plot(r, u_normalized, label=f'E = {E:.2f}')
    plt.plot(r, u_normalized**2, '--', label=f'Probability density E = {E:.2f}')

plt.xlabel('r')
plt.ylabel('u(r) and Probability Density')
plt.legend()
plt.title('Radial Wave Function and Probability Density')
plt.show()

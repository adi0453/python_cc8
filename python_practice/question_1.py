import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Given values
m = 1  # mass
V0 = 20  # potential well height
a = 2  # width of the well
hbar = 1  # 1/hbar^2

# Define the transcendental equation
def lhs(E):
    return np.sqrt(u0**2 - v_values**2)

def rhs(E):
    return v_values * np.tan(v_values)

# Energy range
E = np.linspace(0.1, V0 , 1000)

u0 = np.sqrt(m*(a**2)*V0/2*hbar**2)

def v(E):
    return a/2*np.sqrt(2*m*E/hbar**2)
    
v_values = v(E)


# Plot the functions
plt.plot(E, lhs(E), label=r'$\sqrt{\frac{V_0 - E}{E}}$')
plt.plot(E, rhs(E), label=r'$\tan(\sqrt{\frac{2mE}{\hbar^2}}a)$')
plt.xlim(0, V0)
plt.ylim(-10, 10)
plt.xlabel('Energy (E)')
plt.ylabel('Function values')
plt.title('Graphical Solution of the Transcendental Equation')
plt.legend()
plt.grid(True)
plt.show()

def V(x):
    if x > a/2:
        return V0
    else:
        return 0
x = np.linspace(-a, a, 1000)
Vx = np.array([V(xi) for xi in x])

dx = x[1] - x[0]

def solve_schrodinger_equation(E):
    psi = np.zeros_like(x)
    psi[0] = 1e-5
    psi[1] = 1e-5

    for i in range(1, len(x) - 1):
        psi[i+1] = 2*psi[i] - psi[i-1] + dx**2*(2*m/hbar**2)*(Vx[i] - E)*psi[i]

    def integrand(xi):
        # Interpolate psi at xi
        index = int((xi + a) / dx)
        return psi[index]**2
    
    # Numerical integration to compute the normalization constant
    norm, _ = quad(integrand, -a, a)
    
    # Normalize the wave function
    norm_factor = np.sqrt(norm)
    return psi / norm_factor

energies = [1.5, 4.5]

for E in energies:
    psi = solve_schrodinger_equation(E)
    plt.plot(x, psi)
plt.show()

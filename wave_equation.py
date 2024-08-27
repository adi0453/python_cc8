import numpy as np
import matplotlib.pyplot as plt

# Define grid parameters
nx = 100  # Number of spatial points
nt = 300  # Number of time steps
x = np.linspace(0, 1, nx)
dx = x[1] - x[0]
dt = 0.01  # Time step

# Initialize wave field arrays
U = np.zeros((nx, nt))  # U(x,t)

# Initial condition: U(x,0) = sin(pi*x)
U[:, 0] = np.sin(np.pi * x)

# Apply the wave equation iteratively
for n in range(0, nt-1):
    for i in range(1, nx-1):
        U[i, n+1] = 2*U[i, n] - U[i, n-1] + (dt/dx)**2 * (U[i+1, n] - 2*U[i, n] + U[i-1, n])

# Plot the results at different time steps
plt.plot(x, U[:, 0], label='t=0')
plt.plot(x, U[:, int(nt/3)], label='t=0.5')
plt.plot(x, U[:, int(2*nt/3)], label='t=1.0')
plt.plot(x, U[:, -1], label='t=1.5')

plt.xlabel('x')
plt.ylabel('U(x,t)')
plt.legend()
plt.title('1D Wave Equation')
plt.show()

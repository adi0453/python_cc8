import numpy as np
import matplotlib.pyplot as plt

# Parameters
D = 1.0  # Diffusion coefficient
nx = 100  # Number of spatial points
nt = 200  # Number of time steps
x = np.linspace(0, 10, nx)  # Spatial domain
dx = x[1] - x[0]  # Spatial step size
dt = 0.01  # Time step

# Initialize the temperature field
U = np.zeros((nx, nt))  # U(x,t)

# Apply boundary conditions
U[0, :] = 100  # U(0, t) = 100
U[-1, :] = 100  # U(10, t) = 100

# Time-stepping loop to solve the diffusion equation
for n in range(0, nt-1):
    for i in range(1, nx-1):
        U[i, n+1] = U[i, n] + D * (dt/dx**2) * (U[i+1, n] - 2*U[i, n] + U[i-1, n])

# Plot the results at different time steps
plt.plot(x, U[:, 0], label='t=0')
plt.plot(x, U[:, int(nt/4)], label='t=1.25')
plt.plot(x, U[:, int(nt/2)], label='t=2.5')
plt.plot(x, U[:, -1], label='t=5.0')

plt.xlabel('x')
plt.ylabel('U(x,t)')
plt.legend()
plt.title('1D Diffusion Equation')
plt.show()

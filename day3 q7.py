import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

def exact_sol(omega):
    return F/(np.sqrt((k-omega**2)**2 + gamma**2 * omega**2))

def harmonic_motion(u, t, omega):
    x, y = u
    f1 = y
    f2 = F*np.cos(omega*t)-gamma*y-k*x
    return [f1, f2]

gamma, k, F = 0.1, 2, 1
u0 = [0, 0.001]
t = np.linspace(0, 100, 1000)
omega_values = np.linspace(0.1, 3, 100)
amplitudes = []
for omega in omega_values:
    sol = odeint(harmonic_motion, u0, t, args=(omega, ))
    x_steady_state = sol[:,0][int(0.9 * len(t)):]
    amplitude = np.max(np.abs(x_steady_state))
    amplitudes.append(amplitude)
    
plt.plot(omega_values, amplitudes)
plt.show()

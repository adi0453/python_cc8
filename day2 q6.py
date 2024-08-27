import numpy as np
from scipy.special import legendre as P
import matplotlib.pyplot as plt

#A
x = np.linspace(-1, 1, 10)

def lhs(x, n):
    return (2*n+1)*x*P(n)(x)

def rhs(x, n):
    return (n + 1)*P(n + 1)(x) + n*P(n - 1)(x)
tol = 0.0001
for n in range(1, 5):
    LHS = lhs(x, n)
    RHS = rhs(x, n)
    plt.plot(x, RHS, color="blue")
    plt.plot(x, LHS, color="red")
#plt.show()
    

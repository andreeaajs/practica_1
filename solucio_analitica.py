import numpy as np
import matplotlib.pyplot as plt

t_a = 0.025
T_a = 0.459
N= 101 
x = np.linspace(0, 1, N)

def sumatori(x,t,N = 2000): 
    valor = 0
    for n in range(N):
        valor = valor + (1-np.exp(-((2*n+1)*np.pi)**2 * t)) * np.sin((2*n+1)*np.pi*x) / (2*n+1)**3
    return valor

def T(x,t):
    T = np.zeros(N)
    for i in range(N):
        valol = T_a + (4/(np.pi**3)) * sumatori(x[i],t_a,10**4)
        T[i] = valol
    return T

print(T(x,t_a))
print(x)

plt.figure()
plt.plot(x,T(x,t_a))
plt.show()



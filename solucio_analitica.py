import numpy as np
import matplotlib.pyplot as plt

t_a = 0.025
sigma, V_pic, L = 0.472, 40.0, 0.02
kappa = 0.56
Tc = 36.5

# Calculem la potència externa, la temperatura característica i les temperatures adimensionals.
P_ext = sigma*V_pic**2/(2*L**2)
dTc   = P_ext*L**2/kappa
beta = (Tc + 273.15)/dTc
beta_C = Tc/dTc     #això és necessari?????????

# Discretitzem el domini espacial en N punts equidistants.                      
N= 101
x = np.linspace(0, 1, N)

# Calculem la suma parcial de la sèrie de Fourier.
def sumatori(x,t,N = 2000): 
    valor = 0
    for n in range(N):
        valor = valor + (1-np.exp(-((2*n+1)*np.pi)**2 * t)) * np.sin((2*n+1)*np.pi*x) / (2*n+1)**3
    return valor

# Calculem la temperatura analítica en cada punt espacial.
def T(x,t):
    T = np.zeros(N)
    for i in range(N):
        valol = beta + (4/(np.pi**3)) * sumatori(x[i],t_a,10**4)
        T[i] = valol
    return T

print(T(x, t_a))    #aixo no cal!!!!
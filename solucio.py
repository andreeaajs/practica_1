import numpy as np
import _euler_explicit_
from _euler_explicit_ import euler_explicit

N = 101
x = np.linspace(0,1, N)
dx = x[1] - x[0]
dt = 0.25 * dx**2                  

malalt = (x >= 0.375) & (x <= 0.625)
sa = ~malalt

sigma, V_pic, L = 0.472, 40.0, 0.02      
kappa = 0.56                             
Tc = 36.5                                
P_ext = sigma*V_pic**2/(2*L**2)      
dTc   = P_ext*L**2/kappa                 

#Definim les temperatures màxima (80ºC) i mínima (90ºC) adimensional
T_50 = (50 - Tc)/dTc                     
T_80 = (80 - Tc)/dTc 


def simulacio(t_off, t_final):
#Vector de temperatura
    T = np.zeros(N)
    n_passos = int(t_final/dt)
    n_off = int(t_off/dt)

    temps = 0.0
    temps_sa = 0.0
    temps_Tmax = 0.0

    for n in range(n_passos):
        if n < n_off:
            q = 1.0
        else:
            q = 0.0

        #Euler explícit
        T_nou = T.copy()
        for i in range(1, N-1):
            T_nou[i] = T[i] + _euler_explicit_.euler_explicit(dt, dx, T[i+1], T[i], T[i-1])
        T = np.vstack([T, T_nou])

        n_malalts = 0
        n_malalts_banda = 0 
        sa_supera_50 = False 

        for i in range(N):

            if malalt[i]:
                n_malalts += 1
                if T[i] >= T_50 and T[i] < T_80:
                    n_malalts_banda += 1 
            
            else:
                if T[i] >= T_50:
                    sa_supera_50 += 1 



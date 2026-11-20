import numpy as np

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

#Vector de temperatura
T = np.zeros(N)
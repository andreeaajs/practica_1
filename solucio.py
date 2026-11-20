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
print(T_80, T_50)


vector = np.zeros(N)

def condicions(T):
    T_esquerra = T[0:38] 
    T_malalt = T[38:63]
    T_dreta = T[63:101]

    if np.max(T_esquerra) < T_50 and np.max(T_dreta) < T_50 and np.max(T_malalt) < T_80:
        return True 
    else: 
        return False 

temps = [vector.copy()]


while condicions(vector) == True:
    nova_fila = vector.copy()
    for i in range(1, N-1):
        nova_fila[i] = vector[i] + _euler_explicit_.euler_explicit(dt, dx, vector[i+1], vector[i], vector[i-1])
    vector = nova_fila
    temps.append(vector.copy())
    

temps = np.array(temps)
bo = temps[:-1]
print('temps òptim:', (len(bo) - 1)*dt)


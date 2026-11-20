import numpy as np
import _euler_explicit_
from _euler_explicit_ import euler_explicit

N = 101
x = np.linspace(0,1, N)
dx = x[1] - x[0]
dt = 0.25 * dx**2                  
T_a = 0.459

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


matriu_1 = np.array([[T_a]*N], dtype=float)
vector = matriu_1[-1]

def condicions(Tmax, Tmin):
    T_esquerra = vector[0:38] 
    T_dreta = vector[38:62]
    T_malalt = vector[62:101]

    if np.max(T_esquerra) < Tmax and np.max(T_dreta) < Tmax and np.max(T_malalt) < 80:
        return True 
    else: 
        return False 



while condicions(T_80, T_50) == True:
    nova_fila = matriu_1[i].copy()
    for i in range(1, N-1):
        nova_fila[i] = matriu_1[i] + _euler_explicit_.euler_explicit(dt, dx, matriu_1[i+1], matriu_1[i], matriu_1[i-1])
    matriu_1 = np.vstack([matriu_1, nova_fila])

vector_temp = matriu_1[-1]  

print('hola')
print(matriu_1)






















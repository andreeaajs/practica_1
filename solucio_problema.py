import numpy as np
import matplotlib.pyplot as plt
import _euler_explicit_
from _euler_explicit_ import euler_explicit

#Definim els passos d'espai i temps 
N = 101
x = np.linspace(0,1, N)
dx = x[1] - x[0]
dt = 0.25 * dx**2                  


#Definim les constants que ens donen al guió i calculem Pext
sigma, V_pic, L = 0.472, 40.0, 0.02      
kappa = 0.56                             
Tc = 36.5                                
P_ext = sigma*V_pic**2/(2*L**2)      
dTc   = P_ext*L**2/kappa    
cv, rho = 3698, 1081
alpha = kappa/(cv*rho)

#Definim les temperatures màxima (80ºC) i mínima (90ºC) adimensional (T=0 correspon a 36.6ºC)
T_50 = (50 - Tc)/dTc                     
T_80 = (80 - Tc)/dTc 

#Aquest vector és el que utilitzarem per fer la nostra iteració 
vector = np.zeros(N)

#Definim una funció que ens torni True mentre es compleixin les condicions del problema i False si no es compleixen
def condicions(T):
    T_esquerra = T[0:38] #Diferenciem la part del vector T en la regió de la esquerra de la part malalta, la part malalta i la part de la dreta d'aquesta
    T_malalt = T[38:63]
    T_dreta = T[63:101]

    if np.max(T_esquerra) < T_50 and np.max(T_dreta) < T_50 and np.max(T_malalt) < T_80:
        return True 
    else: 
        return False 

temps = [vector.copy()] #Cream un vector per guardar els vectors de tempertatures a cada pas de temps                                                                                                                                                                                              

#Aqui fem el bucle que calcula 
while condicions(vector) == True:
    nova_fila = vector.copy() #còpia per mantenir les fronteres
    for i in range(1, N-1): 
        #Fem euler explícit
        nova_fila[i] = vector[i] + _euler_explicit_.euler_explicit(dt, dx, vector[i+1], vector[i], vector[i-1])
    vector = nova_fila   # el nou vector passa a ser l'estat actual
    temps.append(vector.copy())
    

temps = np.array(temps) 
fila_temps_optim = temps[:-1] #Obtenim l'ultim pas guardat
temps_optim = (len(fila_temps_optim) - 1)*dt #El temps (adimensional) de l'últim estat vàlid
temps_opt_s = temps_optim * L**2/alpha #Passem el temps adimensional a segons

# Perfil de temperatura a l'instant òptim
T_final_adim = fila_temps_optim[-1]          # última fila de la matriu de temps
T_final = Tc + dTc * T_final_adim            # passem a ºC: T = Tc + T̂·dTc
x_cm = x * L * 100                           # passem x̂ a cm


print('Temps òptim adimensional:', temps_optim)
print('Temps òptim en segons:', temps_opt_s, 's')
print(T_50, T_80)


#Representació de la distribució de temperatures en l'intant que s'apaga la font
plt.figure()
plt.plot(x_cm, T_final)
plt.show()
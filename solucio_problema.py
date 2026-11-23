import numpy as np
import matplotlib.pyplot as plt
import _euler_explicit_


# Configurem l'estil i el format dels gràfics.

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm", 
    "font.size": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 8,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.2,
    "legend.frameon": False,
    "figure.dpi": 120,
    "savefig.bbox": "tight",
})

# Definim la discretització i els passos espacials i temporals.
N = 101
x = np.linspace(0,1, N)
dx = x[1] - x[0]
dt = 0.25 * dx**2

# Definim els paràmetres físics i calculem la potència externa,
# l'increment característic de temperatura i la difusivitat tèrmica.
sigma, V_pic, L = 0.472, 40.0, 0.02      
kappa = 0.56                             
Tc = 36.5                                
P_ext = sigma*V_pic**2/(2*L**2)      
dTc   = P_ext*L**2/kappa    
cv, rho = 3698, 1081
alpha = kappa/(cv*rho)

# Definim els límits de temperatura adimensionals:
# 50ºC per al teixit sa i 80ºC per al teixit malalt
# (T=0 correspon a 36.6ºC).
T_50 = (50 - Tc)/dTc                     
T_80 = (80 - Tc)/dTc 

# Aquest vector és el que utilitzarem per fer la nostra iteració.
vector = np.zeros(N)

#Definim una funció que ens torna True mentre es compleixin les condicions del problema i False si no es compleixen.
def condicions(T):
    T_esquerra = T[0:38] #Diferenciem la part del vector T en la regió de la esquerra de la part malalta, la part malalta i la part de la dreta d'aquesta.
    T_malalt = T[38:63]
    T_dreta = T[63:101]

    if np.max(T_esquerra) < T_50 and np.max(T_dreta) < T_50 and np.max(T_malalt) < T_80:
        return True 
    else: 
        return False 

# Cream un vector per guardar els vectors de tempertatures a cada pas de temps.
temps = [vector.copy()]

# Evolucionem la temperatura mentre es compleixin els límits establerts.
while condicions(vector) == True:
    nova_fila = vector.copy() # Còpia per mantenir les fronteres.
    for i in range(1, N-1): 
        #Fem euler explícit.
        nova_fila[i] = vector[i] + _euler_explicit_.euler_explicit(dt, dx, vector[i+1], vector[i], vector[i-1])
    vector = nova_fila   # El nou vector passa a ser l'estat actual.
    temps.append(vector.copy())
    
# Seleccionem l'últim estat que compleix les condicions
# i calculem el temps òptim en unitats adimensionals i en segons.
temps = np.array(temps) 
fila_temps_optim = temps[:-1] # Obtenim l'ultim pas guardat.
temps_optim = (len(fila_temps_optim) - 1)*dt # El temps (adimensional) de l'últim estat vàlid.
temps_opt_s = temps_optim * L**2/alpha # Passem el temps adimensional a segons.

# Perfil de temperatura a l'instant òptim.
T_final_adim = fila_temps_optim[-1]          # Última fila de la matriu de temps.
T_final = Tc + dTc * T_final_adim            # Passem a ºC: T = Tc + T̂·dTc
x_cm = x * L * 100                           # Passem x̂ a cm


print('Temps òptim adimensional:', temps_optim)
print('Temps òptim en segons:', temps_opt_s, 's')
print('Temperatura de la forntera del teixit malalt i el sa:')
print('Frontera teixit sa:', T_final[37], 'ºC.', 'Frontera teixit malalt', T_final[38], 'ºC.')
print('Temperatura al mig del teixit malalt:', T_final[50]) 
# La temperatura del teixit malalt només supera lleugerament els 50 ºC, 
# ja que la temperatura a la frontera amb el teixit sa és molt pròxima a aquest límit.

#Representació de la distribució de temperatures en l'instant que s'apaga la font.
plt.figure()
plt.plot(x_cm, T_final)
plt.axvspan(0.375*L*100, 0.625*L*100, color='tab:blue', alpha=0.2, label='Teixit malalt')
plt.axhline(50, color = 'orange', linestyle = '--', label = '50ºC')
plt.xlabel('Distància x (cm)')
plt.ylabel('Temperatura (ºC)')
plt.legend()
plt.show()
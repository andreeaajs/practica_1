import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs
from solucio_analitica import T as T_analitica

# Perque els grafics quedin millor

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



#PARÀMETRES
N = 101
x = np.linspace(0, 1, N)
dx = x[1] - x[0]
t_a = 0.025

#aquests calculen la T_a (normalitzada) correctament
sigma, V_pic, L = 0.472, 40.0, 0.02      
kappa = 0.56                             
Tc = 36.5                                
P_ext = sigma*V_pic**2/(2*L**2)        
dTc   = P_ext*L**2/kappa                 

T_a = (Tc + 273.15)/dTc               


dt_4 = dx**2
dt_5 = 0.5*dx**2
gamma_1 = 1
gamma_2 = 0.5

passos_4 = round(0.025/dt_4)
passos_5 = round(0.025/dt_5)

# Començarem amb una matriu plena amb T_a i anirem reescribint les files.
Matriu_T_4 = np.full((passos_4+1, N), T_a) # La matriu ha de fer passos_4+1 passos de temps ja que la primera fila és el temps inicial.



# Definim els parametres off i diag i un vector que guardarà el nombre d'iteracions
off_4 = gamma_1/2
diag_4 = 1 + gamma_1
iters_4 = [0]*(passos_4+1)
b_1 = np.zeros(N)

# Comencem el bucle temporal

for t in range(passos_4):
    # Definim el vector de termes independents b que va canviant a cada pas
    for j in range(1, N-1):
        b_1[j] = Matriu_T_4[t][j] + (gamma_1/2)*(Matriu_T_4[t][j+1] - 2*Matriu_T_4[t][j] + Matriu_T_4[t][j-1]) + dt_4

    T_p = Matriu_T_4[t].copy()  # punt de partida per a Gauss-Seidel
    nou, its = gs.gauss_seidel(T_p, b_1, off_4, diag_4)
    Matriu_T_4[t+1] = nou
    iters_4[t+1] = its


vector_4 = Matriu_T_4[-1]

error_4 = np.abs(vector_4 - T_analitica(x,t_a))


#print(vector_4)

plt.figure()
plt.scatter(x, vector_4, color = '#1B4F72', s=7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')   
plt.title('Distribució de temperatura amb Crank-Nicolson per $\gamma=1$')

plt.figure()
plt.plot(x,error_4, color = '#1B4F72')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')  
plt.title('Error numèric amb Crank-Nicolson $\gamma=1$')


# Repetim tot l'anterior per a dt_5 ------------------------------------------------------

# Començarem amb una matriu plena amb T_a i anirem reescribint les files.
Matriu_T_5 = np.full((passos_5+1, N), T_a) # La matriu ha de fer passos_4+1 passos de temps ja que la primera fila és el temps inicial.



# Definim els parametres off i diag i un vector que guardarà el nombre d'iteracions
off_5 = gamma_2/2
diag_5 = 1 + gamma_2
iters_5 = [0]*(passos_5+1)
b_2 = np.zeros(N)

# Comencem el bucle temporal

for t in range(passos_5):
    # Definim el vector de termes independents b que va canviant a cada pas
    for j in range(1, N-1):
        b_2[j] = Matriu_T_5[t][j] + (gamma_2/2)*(Matriu_T_5[t][j+1] - 2*Matriu_T_5[t][j] + Matriu_T_5[t][j-1]) + dt_5

    T_p = Matriu_T_5[t].copy()  # punt de partida per a Gauss-Seidel
    nou, its = gs.gauss_seidel(T_p, b_2, off_5, diag_5)
    Matriu_T_5[t+1] = nou
    iters_5[t+1] = its


vector_5 = Matriu_T_5[-1]
print(vector_5)

plt.figure()
plt.scatter(x, vector_5, color = '#922B21', s=7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')      
plt.title('Distribució de temperatura amb Crank-Nicolson $\gamma=0.5$')

err_5 = np.abs(T_analitica(x,t_a)-vector_5)

plt.figure()
plt.plot(x,err_5, color = '#922B21')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')   
plt.title('Error numèric amb Crank-Nicolson $\gamma=0.5$')
plt.show()


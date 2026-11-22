import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs
from Solucio_analitica import T as T_analitica
from Solucio_problema import beta


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
T_a = beta
t_a = 0.025
x = np.linspace(0, 1, N)
dx = x[1] - x[0]

#-----------------------------EULER IMPLÍCIT---------------------------------

# Els paràmetres d'aquest mètode són:
dt_4 = dx**2
dt_5 = 0.5*dx**2
gamma_1 = 1
gamma_2 = 0.5

passos_4 = round(0.025/dt_4)
passos_5 = round(0.025/dt_5)

# Començarem amb una matriu plena amb T_a i anirem reescribint les files.

Matriu_T_4 = np.full((passos_4+1, N), T_a) # La matriu ha de fer passos_4+1 passos de temps 
# ja que la primera fila és el temps inicial.

Matriu_T_5 = np.full((passos_5+1, N), T_a)





# Definim els parametres off i diag i un vector que guardarà el nombre d'iteracions
off_4 = gamma_1
diag_4 = 1 + 2*gamma_1
iters_4 = [0]*(passos_4+1)

off_5 = gamma_2
diag_5 = 1 + 2*gamma_2
iters_5 = [0]*(passos_5+1)

# Comencem el bucle temporal

for t in range(passos_4):
    # Definim el vector de termes independents b que va canviant a cada pas
    b = Matriu_T_4[t] + dt_4

    T = Matriu_T_4[t].copy()  # punt de partida per a Gauss-Seidel
    nou, its = gs.gauss_seidel(T, b, off_4, diag_4)
    Matriu_T_4[t+1] = nou
    iters_4[t+1] = its

# Repetim tot l'anterior per a dt_5 ------------------------------------------------------


for t in range(passos_5):
    # Definim el vector de termes independents b que va canviant a cada pas
    b = Matriu_T_5[t] + dt_5

    T = Matriu_T_5[t].copy()  # punt de partida per a Gauss-Seidel
    nou, its = gs.gauss_seidel(T, b, off_5, diag_5)
    Matriu_T_5[t+1] = nou
    iters_5[t+1] = its

# Representació gràfica dels resultats

vector_4 = Matriu_T_4[-1]

error_4 = np.abs(vector_4-T_analitica(x,t_a))

plt.figure()
plt.scatter(x, vector_4, color='#1B4F72', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}$')   
plt.title('Distribució de temperatura a $\hat{t}=0.025$ amb Euler implícit per $\Delta \hat{t}=\Delta \hat{x}^2$')

plt.figure()
plt.plot(x, error_4, color='#1B4F72')
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}$')   
plt.title('Error numèric per $\Delta \hat{t}=\Delta \hat{x}^2$')




vector_5 = Matriu_T_5[-1]

error_5 = np.abs(vector_5 - T_analitica(x,t_a))

plt.figure()
plt.scatter(x, vector_5, color='#922B21', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}$')   
plt.title('Distribució de temperatura a $\hat{t}=0.025$ amb Euler implícit per $\Delta \hat{t}=0.5\Delta \hat{x}^2$')

plt.figure()
plt.plot(x, error_5, color='#922B21')
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}$')   
plt.title('Error numèric per $\Delta \hat{t}=0.5\Delta \hat{x}^2$')
plt.show()





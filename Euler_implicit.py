import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs
from Solucio_analitica import T , t_a

#PARÀMETRES

N = 101
T_a = 0.459
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

plt.scatter(x, vector_4, label='dt = dx^2', color='green')
plt.xlabel('x')
plt.ylabel('Temperatura')   
plt.title('Distribució de temperatura a t=0.025 amb Euler implícit')
plt.legend()
plt.show()


vector_5 = Matriu_T_5[-1]

plt.scatter(x, vector_5, label='dt = 0.5*dx^2', color='blue')
plt.xlabel('x')
plt.ylabel('Temperatura')   
plt.title('Distribució de temperatura a t=0.025 amb Euler implícit')
plt.legend()
plt.show()

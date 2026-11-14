import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs

#PARÀMETRES

N = 101
T_a = 0.459
x = np.linspace(0, 1, N)
dx = x[1] - x[0]

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

    T = Matriu_T_4[t].copy()  # punt de partida per a Gauss-Seidel
    nou, its = gs.gauss_seidel(T, b_1, off_4, diag_4)
    Matriu_T_4[t+1] = nou
    iters_4[t+1] = its


vector_4 = Matriu_T_4[-1]
print(vector_4)

#plt.scatter(x, vector_4, label='dt = dx^2', color='green')
plt.xlabel('x')
plt.ylabel('Temperatura')   
plt.title('Distribució de temperatura a t=0.025 amb Euler implícit')
plt.legend()
#plt.show()

# Repetim tot l'anterior per a dt_5 ------------------------------------------------------

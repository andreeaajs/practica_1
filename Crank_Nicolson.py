import numpy as np
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs
from solucio_analitica import T as T_analitica

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


#PARÀMETRES:

N = 101                     # Nombre de punts de la malla espacial.
x = np.linspace(0, 1, N)    # Discretització espacial.
dx = x[1] - x[0]            # Pas espacial.
t_a = 0.025                 # Temps final.

# Calculem la temperatura inicial adimensional.
sigma, V_pic, L = 0.472, 40.0, 0.02
kappa = 0.56
Tc = 36.5
P_ext = sigma*V_pic**2/(2*L**2)
dTc   = P_ext*L**2/kappa
T_a = (Tc + 273.15)/dTc

# Definim els passos temporals i els paràmetres dels dos casos.
dt_4 = dx**2
dt_5 = 0.5*dx**2
gamma_1 = 1
gamma_2 = 0.5

# Nombre de passos temporals.
passos_4 = round(t_a/dt_4)
passos_5 = round(t_a/dt_5)


# CRANCK NICOLSON:

# Comencem amb una matriu plena amb T_a i anirem reescribint les files.
# Afegim una fila per incloure l'instant inicial.
Matriu_T_4 = np.full((passos_4+1, N), T_a)

# Definim els parametres d eles diagonals (off i diag).
off_4 = gamma_1/2
diag_4 = 1 + gamma_1

# Vector que emmagatzemarà les iteracions.
iters_4 = [0]*(passos_4+1)

# Vector de termes independents.
b_1 = np.zeros(N)

# Evolució temporal.
for t in range(passos_4):
    for j in range(1, N-1):     # Definim el vector de termes independents b que va canviant a cada pas.
        b_1[j] = Matriu_T_4[t][j] + (gamma_1/2)*(Matriu_T_4[t][j+1] - 2*Matriu_T_4[t][j] + Matriu_T_4[t][j-1]) + dt_4
    T_p = Matriu_T_4[t].copy()  # Punt de partida per a Gauss-Seidel.
    nou, its = gs.gauss_seidel(T_p, b_1, off_4, diag_4)     # Resolem el sistema lineal amb Gauss-Seidel.
    # Emmagatzemem la nova temperatura i les iteracions.
    Matriu_T_4[t+1] = nou
    iters_4[t+1] = its

# Distribució final i error absolut.
vector_4 = Matriu_T_4[-1]
error_4 = np.abs(vector_4 - T_analitica(x,t_a))


# Representem la distribució de temperatura i l'error del primer cas.
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


# Repetim tot l'anterior per a dt_5.
# Comencem amb una matriu plena amb T_a i anirem reescribint les files.
# Afegim una fila per incloure l'instant inicial.
Matriu_T_5 = np.full((passos_5+1, N), T_a)

# Definim els parametres d eles diagonals (off i diag).
off_5 = gamma_2/2
diag_5 = 1 + gamma_2

# Vector que emmagatzemarà les iteracions.
iters_5 = [0]*(passos_5+1)

# Vector de termes independents.
b_2 = np.zeros(N)

# Evolució temporal.
for t in range(passos_5):
    for j in range(1, N-1):         # Definim el vector de termes independents b que va canviant a cada pas.
        b_2[j] = Matriu_T_5[t][j] + (gamma_2/2)*(Matriu_T_5[t][j+1] - 2*Matriu_T_5[t][j] + Matriu_T_5[t][j-1]) + dt_5
    T_p = Matriu_T_5[t].copy()  # Punt de partida per a Gauss-Seidel.
    nou, its = gs.gauss_seidel(T_p, b_2, off_5, diag_5)     # Resolem el sistema lineal amb Gauss-Seidel.
    # Emmagatzemem la nova temperatura i les iteracions.
    Matriu_T_5[t+1] = nou
    iters_5[t+1] = its

# Distribució final i error absolut.
vector_5 = Matriu_T_5[-1]
err_5 = np.abs(T_analitica(x,t_a)-vector_5)

print(vector_5)

# Representem la distribució de temperatura i l'error del segon cas.
plt.figure()
plt.scatter(x, vector_5, color = '#922B21', s=7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')      
plt.title('Distribució de temperatura amb Crank-Nicolson $\gamma=0.5$')

plt.figure()
plt.plot(x,err_5, color = '#922B21')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')   
plt.title('Error numèric amb Crank-Nicolson $\gamma=0.5$')

plt.show()
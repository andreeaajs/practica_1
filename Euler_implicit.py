import numpy as np
import matplotlib.pyplot as plt
import _gauss_seidel_ as gs
from solucio_analitica import T as T_analitica
from solucio_analitica import beta


# Configurem l'estil i el format dels gràfics.

plt.rcParams.update({
    "font.family": "serif", # font de la lletra
    "mathtext.fontset": "cm", # les formules es fan amb "computer modern" igual que LaTeX
    "font.size": 10, # mida de la lletra
    "axes.labelsize": 10, # mida de les etiquetes
    "legend.fontsize": 8, # mida del text de la llegenda
    "xtick.labelsize": 9, # mida dels números de les marques dels eixos
    "ytick.labelsize": 9, 
    "xtick.direction": "in", "ytick.direction": "in", # Les marques cap a "dins"
    "xtick.top": True, "ytick.right": True, # dibuxa també marques als costats dret i superior
    "xtick.minor.visible": True, "ytick.minor.visible": True, # Afegeix marques menors
    "axes.linewidth": 0.8, # Gruix de la caixa dels eixos
    "lines.linewidth": 1.2, # gruix de les corbes de plot
    "legend.frameon": False, # Treu la quadricula
    "figure.dpi": 120, # marca la resolució
    "savefig.bbox": "tight", # si es guarda la figura, retalla l'espai en blanc
})

#PARÀMETRES:

N = 101                     # Nombre de punts de la malla espacial.
T_a = beta                  # Temperatura inicial.
t_a = 0.025                 # Temps final.
x = np.linspace(0, 1, N)    # Discretització espacial.
dx = x[1] - x[0]            # Pas espacial.


# EULER IMPLÍCIT:

# Definim els passos temporals i els paràmetres dels dos casos.
dt_4 = dx**2
dt_5 = 0.5*dx**2
gamma_1 = 1
gamma_2 = 0.5

# Calculem el nombre de passos temporals de cada cas.
passos_4 = round(0.025/dt_4)
passos_5 = round(0.025/dt_5)

# Definim les matrius, que comencen amb la temperatura inicial.
# Afegim una fila al nombre de passos per incloure l'instant inicial.
Matriu_T_4 = np.full((passos_4+1, N), T_a)
Matriu_T_5 = np.full((passos_5+1, N), T_a)

# Definim els coeficients de les diagonals de cada sistema.
off_4 = gamma_1
diag_4 = 1 + 2*gamma_1
off_5 = gamma_2
diag_5 = 1 + 2*gamma_2

# Definim els vectors que emmagatzemaran les iteracions de Gauss-Seidel.
iters_4 = [0]*(passos_4+1)
iters_5 = [0]*(passos_5+1)

# Evolució temporal.
for t in range(passos_4):
    b = Matriu_T_4[t] + dt_4    # Definim el vector de termes independents b que va canviant a cada pas.
    T = Matriu_T_4[t].copy()  # Punt de partida per a Gauss-Seidel.
    nou, its = gs.gauss_seidel(T, b, off_4, diag_4)     # Resolem el sistema lineal amb Gauss-Seidel.
    # Guardem la nova temperatura i les iteracions.
    Matriu_T_4[t+1] = nou
    iters_4[t+1] = its

# Repetim tot l'anterior per a dt_5.
for t in range(passos_5):
    b = Matriu_T_5[t] + dt_5
    T = Matriu_T_5[t].copy()
    nou, its = gs.gauss_seidel(T, b, off_5, diag_5)
    Matriu_T_5[t+1] = nou
    iters_5[t+1] = its


# Representació gràfica dels resultats:

# Distribució final i error absolut.
vector_4 = Matriu_T_4[-1]
error_4 = np.abs(T_analitica(x, (passos_4*dt_4)) - vector_4)

# Representem la distribució de temperatura i l'error del primer cas.
plt.figure()
plt.scatter(x, vector_4, color='#1B4F72', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')   

plt.figure()
plt.plot(x, error_4, color='#1B4F72')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')

# Distribució final i error absolut.
vector_5 = Matriu_T_5[-1]
error_5 = np.abs(T_analitica(x, (passos_5*dt_5)) - vector_5)

# Representem la distribució de temperatura i l'error del segon cas.
plt.figure()
plt.scatter(x, vector_5, color='#922B21', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')   

plt.figure()
plt.plot(x, error_5, color='#922B21')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')

plt.show()
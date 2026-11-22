import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
from solucio_analitica import T
from solucio_analitica import beta 
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
t_a = 0.025
T_a = beta
x = np.linspace(0, 1, N)
dx = x[1] - x[0]
dt_1 = 0.51*dx**2
dt_2 = 0.49*dx**2
dt_3 = 0.25*dx**2
gamma_1 = 1
gamma_2 = 0.5

#EULER EXPLICIT

passos_t1 = round(t_a/dt_1)
passos_t2 = round(t_a/dt_2)
passos_t3 = round(t_a/dt_3)


matriu_1 = np.array([[T_a]*N], dtype=float)

for t in range(passos_t1):
    nova_fila = matriu_1[t].copy()
    for i in range(1, N-1):
        nova_fila[i] = matriu_1[t][i] + _euler_explicit_.euler_explicit(dt_1, dx, matriu_1[t][i+1], matriu_1[t][i], matriu_1[t][i-1])
    matriu_1 = np.vstack([matriu_1, nova_fila])


matriu_2 = np.array([[T_a]*N], dtype=float)

for t in range(passos_t2):
    nova_fila = matriu_2[t].copy()
    for i in range(1, N-1):
        nova_fila[i] = matriu_2[t][i] + _euler_explicit_.euler_explicit(dt_2, dx, matriu_2[t][i+1], matriu_2[t][i], matriu_2[t][i-1])
    matriu_2 = np.vstack([matriu_2, nova_fila])

matriu_3 = np.array([[T_a]*N], dtype=float)

for t in range(passos_t3):
    nova_fila = matriu_3[t].copy()
    for i in range(1, N-1):
        nova_fila[i] = matriu_3[t][i] + _euler_explicit_.euler_explicit(dt_3, dx, matriu_3[t][i+1], matriu_3[t][i], matriu_3[t][i-1])
    matriu_3 = np.vstack([matriu_3, nova_fila])


vector_1 = matriu_1[-1]
vector_2 = matriu_2[-1]
vector_3 = matriu_3[-1]

err_1 = np.abs(T(x,t_a) - vector_1)
err_2 = T(x,t_a) - vector_2
err_3 = (T(x,t_a) - vector_3)

# cas 1 i el seu error
plt.figure()
plt.scatter(x, vector_1, color='#1B4F72', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')
plt.title('Distribució de temperatura a $\hat{t}=0.025$ per $\Delta \hat{t} = 0.51\Delta \hat{x}^2$')

plt.figure()
plt.plot(x,err_1, color = '#1B4F72')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')
plt.title('Error numèric per $\Delta \hat{t} = 0.51\Delta \hat{x}^2$')

#cas 2 i el seu error
plt.figure()
plt.scatter(x, vector_2, color='#922B21', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')
plt.title('Distribució de temperatura a $\hat{t}=0.025$ per $\Delta \hat{t} = 0.49\Delta \hat{x}^2$')

plt.figure()
plt.plot(x,err_2, color = '#922B21')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')
plt.title('Error numèric per $\Delta \hat{t} = 0.49\Delta \hat{x}^2$')

#cas 3 i el seu error
plt.figure()
plt.scatter(x, vector_3, color='#1E8449', s = 7)
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')
plt.title('Distribució de temperatura a $\hat{t}=0.025$ per $\Delta \hat{t} = 0.25\Delta \hat{x}^2$')

plt.figure()
plt.plot(x,err_3, color = '#1E8449')
plt.xlabel('$\hat{x}$')
plt.ylabel('Error numèric')
plt.title('Error numèric per $\Delta \hat{t} = 0.25\Delta \hat{x}^2$')

#analítica
plt.figure()
plt.title('Distribució de temperatura a $\hat{t}=0.025$ analítica')
plt.plot(x,T(x,t_a))
plt.xlabel('$\hat{x}$')
plt.ylabel('$\hat{T}(\hat{x})$')

plt.show()
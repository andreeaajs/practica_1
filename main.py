import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt
from solucio_analitica import T


#PARÀMETRES

N = 101
t_a = 0.025
T_a = 0.459
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

plt.figure()
plt.scatter(x, vector_1, label='dt = 0.51*dx^2', color='blue')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('Distribució de temperatura a t=0.025 (adimensional)')
plt.legend()

plt.figure()
plt.scatter(x, vector_2, label='dt = 0.49*dx^2', color='orange')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('Distribució de temperatura a t=0.025 (adimensional)')
plt.legend()

plt.figure()
plt.scatter(x, vector_3, label='dt = 0.25*dx^2', color='green')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('Distribució de temperatura a t=0.025 (adimensional)')
plt.legend()

plt.figure()
plt.title('Distribució de temperatura a t=0.025 (adimensional)')
plt.plot(x,T(x,t_a), label='solució analítica')
plt.legend()

plt.show()
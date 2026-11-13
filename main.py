import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt

N = 101
T_a = 0.459
x = np.linspace(0, 1, N)
dx = x[1] - x[0]
dt_1 = 0.51*dx**2
dt_2 = 0.49*dx**2
dt_3 = 0.25*dx**2



passos_t1 = round(0.025/dt_1)
passos_t2 = round(0.025/dt_2)
passos_t3 = round(0.025/dt_3)



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
plt.title('Temperature distribution at t = 0.025')

plt.figure()
plt.scatter(x, vector_2, label='dt = 0.49*dx^2', color='orange')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('Temperature distribution at t = 0.025')

plt.figure()
plt.scatter(x, vector_3, label='dt = 0.25*dx^2', color='green')
plt.xlabel('x')
plt.ylabel('T(x)')
plt.title('Temperature distribution at t = 0.025')
plt.legend()
plt.show()











import numpy as np
import _euler_explicit_ 
import matplotlib.pyplot as plt

#PARÀMETRES

N = 101
T_a = 0.459
x = np.linspace(0, 1, N)
dx = x[1] - x[0]
dt_4 = dx**2
dt_5 = 0.5*dx**2
gamma_1 = 1
gamma_2 = 0.5

#EULER IMPLÍCIT

passos_4 = round(0.025/dt_4)
passos_5 = round(0.025/dt_5)


T = np.zeros((passos_4+1, N))


def gauss_seidel(T, dt, off, diag, passos):
    T[0] = [T_a]*N
    for t in range(passos):
        err = 0.0
        if t==0:
            b=T[0]+dt
        else:
            b=T[t-1]+dt

        for i in range(1, N-1):
            valor = (b[i] + off*(T[t+1][i-1] + T[t+1][i+1]))/diag 
            err = max(err, abs(valor - T[t+1][i]))
            T[t+1][i]=valor
    return T

#Euler implicit

off = gamma_1 
diag = 1+2*gamma_1

T_resultat = gauss_seidel(T, dt_4, off,diag,passos_4)
print(T_resultat)

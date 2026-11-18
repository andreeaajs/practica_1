# Definim les operacions que fa l'equació d'Euler explícit per calcular cada
# valor de temperatura.

# dt és el pas temporal, dx és el pas espacial, f1 és la temperatura a la posició i+1,
# f0 és la temperatura a la posició i, f_1 és la temperatura a la posició i-1. 

# La funció retorna el valor de la temperatura a la posició i al pas de temps següent.

def euler_explicit(dt, dx, f1, f0, f_1):
    return (dt/dx**2) * (f1 - 2*f0 + f_1) + dt 

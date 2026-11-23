# Definim una funció que calcula l'increment de temperatura
# en un punt de la malla espacial entre dos instants 
# consecutius mitjançant el mètode d'Euler explícit.

# Els paràmetres de la funció són els següents:
#   dt = pas temporal (increment de temps adimensional).
#   dx = pas espacial (increment de posició adimensional).
#   f1 = temperatura al punt espacial següent (j+1).
#   f0 = temperatura al punt espacial actual (j).
#   f_1 = temperatura al punt espacial anterior (j-1).

# La funció retorna l'increment de temperatura al punt j entre els instants i i i+1.

def euler_explicit(dt, dx, f1, f0, f_1):
    return (dt/dx**2) * (f1 - 2*f0 + f_1) + dt
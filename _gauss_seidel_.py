# Definim la funció gauss_seidel per resoldre el sistema d'equacions lineals utilitzant 
# el mètode de Gauss-Seidel.

# La funció té com a paràmetres d'entrada un vector T (que serà la fila de la matriu de 
# temperatures), b que és el vector de termes independents, 
# off que és el valor de la diagonal inferior i superior, 
# diag que és el valor de la diagonal principal, tol que és la tolerància 
# i maxit que és el nombre màxim d'iteracions.

def gauss_seidel(T, b, off, diag, tol=1e-12, maxit=10000):
    for it in range(maxit):
        err = 0.0
        for j in range(1, len(T)-1):
            nou = (b[j] + off*(T[j-1] + T[j+1]))/diag
            err = max(err, abs(nou - T[j]))
            T[j] = nou                     
        if err < tol:
            break
    return T, it

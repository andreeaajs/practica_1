# Definim una funció per resoldre un sistema lineal mitjançant Gauss-Seidel.
# T és el vector de temperatures, b el vector de termes independents,
# off i diag són els elements de les diagonals i tol i maxit
# estableixen la tolerància i el nombre màxim d'iteracions respectivament.

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

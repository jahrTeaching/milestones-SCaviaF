from numpy import linspace, size, zeros, log10
from numpy.linalg import norm
from ODEs.CauchyCode import Cauchy

"Ratio de convergencia"

def Conv(F, t, U0, Scheme,p):

    #p es el numero de puntos de la gráfica que se van a calcular

    #Calculo la primera solución con N divisiones
    N = size(t)
    tf = t[N-1]
    U = Cauchy(F, t, U0, Scheme)

    #Inicializo los vectores
    E = zeros(p)
    Elog = zeros(p)
    Nlog = zeros(p)
    N = 2*N


    for i in range(0,p):

        #Calculo la solución con paso mitad
        t2 = linspace(0, tf, (2**i)*N)
        U2N = Cauchy(F, t2, U0, Scheme)

        #Calculo el ratio de convergencia y su logaritmo
        E[i] = norm((U2N[int((2**i)*N-1),:] - U[int((2**i)*N/2-1),:]))
        Elog[i] = log10(E[i])
        Nlog[i] = log10((2**i)*N)

        U = U2N

    return [Elog, Nlog]

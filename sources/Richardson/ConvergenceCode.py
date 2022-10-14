from numpy import linspace, size, zeros, log10
from numpy.linalg import norm
from ODEs.CauchyCode import Cauchy

def Conv(F, t, U0, Scheme):

    N = size(t)
    tf = t[N-1]
    U = Cauchy(F, t, U0, Scheme)

    p = 8 # Puntos en la gráfica
    E = zeros(p)
    Elog = zeros(p)
    Nlog = zeros(p)
    N = 2*N

    for i in range(0,p):

        t2 = linspace(0, tf, (2**i)*N)
        U2N = Cauchy(F, t2, U0, Scheme)

        E[i] = norm((U2N[int((2**i)*N-1),:] - U[int((2**i)*N/2-1),:]))
        Elog[i] = log10(E[i])
        Nlog[i] = log10((2**i)*N)

        U = U2N
        print(i)

    return [Elog, Nlog]

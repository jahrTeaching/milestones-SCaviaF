from numpy import size, zeros, linspace, log10
from numpy.linalg import norm
from ODEs.CauchyCode import Cauchy

def Conv(F, t, U0, Scheme):

    N = size(t)
    t1 = t
    tf = t1[N-1]
    U = Cauchy(F, t1, U0, Scheme)

    p = 8 # Puntos en la gráfica
    E = zeros(p)
    Elog = zeros(p)
    Nlog = zeros(p)

    for i in range(0,p):

        N = 2*N
        t2 = linspace(0, tf, N)
        U2N = Cauchy(F, t2, U0, Scheme)

        E[i] = norm((U2N[int(N-1),:] - U[int(N/2-1),:]))
        Elog[i] = log10(E[i])
        Nlog[i] = log10(N)

        t1 = t2
        U = U2N
        print(i)

    return [Elog, Nlog]
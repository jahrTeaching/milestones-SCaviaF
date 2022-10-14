from numpy import zeros


"Cauchy Problem"

def Cauchy(F, t, U0, Scheme):

    U = zeros((len(t),len(U0)))
    U[0,:] = U0
    for n in range(len(t)-1):
        U[n+1,:] = Scheme(U[n, :], t[n+1] - t[n], t[n],  F)

    return U

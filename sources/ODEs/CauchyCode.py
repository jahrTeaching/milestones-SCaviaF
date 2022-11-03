from numpy import zeros
import ODEs.SchemesCode as Sch

"Cauchy Problem"

def Cauchy(F, t, U0, Scheme):

    U = zeros((len(t),len(U0)))
    U[0,:] = U0
    dt = t[2] - t[1]

    if Scheme == Sch.LF:
        U[1,:] = U[0,:] + dt*F(U[0,:], t[0])

        for n in range(1,len(t)-1):
            U1 = U[n-1, :]
            U2 = U[n, :]
            U[n+1, :] = Scheme(U2, U1, dt, t[n], F)

    else:
        for n in range(len(t)-1):
            U[n+1,:] = Scheme(U[n, :], t[n+1] - t[n], t[n],  F)

    return U

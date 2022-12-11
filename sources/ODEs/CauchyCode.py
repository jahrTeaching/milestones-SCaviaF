from numpy import zeros
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv, LF

"Cauchy Problem"

def Cauchy(F, t, U0, Scheme):

    U = zeros((len(t),len(U0)))
    U[0,:] = U0
    dt = t[2] - t[1]

    if Scheme == LF:
        U[1,:] = U[0,:] + dt*F(U[0,:], t[0])

        for n in range(1,len(t)-1):
            U1 = U[n-1, :]
            U2 = U[n, :]
            U[n+1, :] = Scheme(U2, U1, dt, t[n], F)

    else:
        for n in range(len(t)-1):
            U[n+1,:] = Scheme(U[n, :], t[n+1] - t[n], t[n],  F)
            if n == 20000:
                print(n)
            if n == 40000:
                print(n)
            if n == 60000:
                print(n)
            if n == 80000:
                print(n)

    return U




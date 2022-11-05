from numpy import linspace, size, zeros
from ODEs.CauchyCode import Cauchy
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv

"Error según Richardson"

def Error(F, t, U0, Scheme):

    N = size(t)
    E = zeros([N,size(U0)])

    t1 = t
    t2 = linspace(0, t[N-1], N*2)

    #Calculo soluciones con dos pasos distintos
    U2 = Cauchy(F, t2, U0, Scheme)
    U1 = Cauchy(F, t1, U0, Scheme)

    #Defino el orden del esquema
    if Scheme == RK4:
        q = 4
    elif Scheme == CN:
        q = 2
    else:
        q = 1

    #Calculo el error
    for i in range(0,N):
        E[i,:] = (U2[2*i,:] - U1[i,:]) / (1 - 1 / (2**q))

    return E
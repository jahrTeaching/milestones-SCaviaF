from numpy import  zeros, array, linspace
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

"Schemes"

def Euler(U, dt, t, F): 

    return U + dt*F(U, t)

def RK4(U, dt, t, F):

    k1 = F(U, t)
    Unew1 = U[:] + dt*k1[:]/2
    k2 = F(Unew1, t)
    Unew2 = U[:] + dt*k2[:]/2
    k3 = F(Unew2, t)
    Unew3 = U[:] + dt*k3[:]
    k4 = F(Unew3, t)

    U = U + dt*(k1 + 2*k2 + 2*k3 + k4)/6

    return U

def CN(U, dt, t, F):

    nextstep = F(U, t)
    def func_CN(x):
        return [x[0] - U[0] - (x[2] + nextstep[0])*dt/2,
                x[1] - U[1] - (x[3] + nextstep[1])*dt/2,
                x[2] - U[2] - (-x[0]/((x[0]**2+x[1]**2)**(3/2)) + nextstep[2])*dt/2,
                x[3] - U[3] - (-x[1]/((x[0]**2+x[1]**2)**(3/2)) + nextstep[3])*dt/2]
    U = fsolve(func_CN, [U])
    "Kepler metido dentro de la función, por lo que solo vale para resolver Kepler y no cualquier función"
    return U

def EulerInv(U, dt, t, F):

    nextstep = F(U, t)
    def func_EI(x):
        return [x[0] - U[0] - (x[2] + nextstep[0])*dt,
                x[1] - U[1] - (x[3] + nextstep[1])*dt,
                x[2] - U[2] - (-x[0]/((x[0]**2+x[1]**2)**(3/2)) + nextstep[2])*dt,
                x[3] - U[3] - (-x[1]/((x[0]**2+x[1]**2)**(3/2)) + nextstep[3])*dt]
    U = fsolve(func_EI, [U])
    "Kepler metido dentro de la función, por lo que solo vale para resolver Kepler y no cualquier función"
    return U



"Function"

def Kepler(U, t):

    return  array([U[2], U[3], -U[0]/((U[0]**2 + U[1]**2)**1.5), -U[1]/((U[0]**2 + U[1]**2)**1.5) ] ) 



"Cauchy Problem"

def Cauchy(F,t,U0,Scheme):

    U = zeros((len(t),len(U0)))
    U[0,:] = U0
    for n in range(len(t)-1):
        U[n+1,:] = Scheme(U[n, :], t[n+1] - t[n], t[n],  F)

    return U



"Run"

tf = 20
N = 200
t = linspace(0, tf, N)
U0 = array([1, 0, 0, 1])

U =  Cauchy(Kepler, t, U0, Euler) 

plt.plot(U[:,0] , U[:,1])
plt.show()

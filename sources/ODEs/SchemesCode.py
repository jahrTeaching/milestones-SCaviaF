from scipy.optimize import fsolve

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

    def func_CN(x):
        return x - U -(F(x,t+dt) + F(U, t))*dt/2
    U = fsolve(func_CN, [U])

    return U

def EulerInv(U, dt, t, F):

    def func_EI(x):

        return x - U - F(x,t+dt)*dt
    U = fsolve(func_EI, [U])

    return U

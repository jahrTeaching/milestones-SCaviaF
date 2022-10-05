from numpy import  zeros, array, linspace
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv
from Functions.KeplerCode import Kepler
from ODEs.CauchyCode import Cauchy


"Run"

tf = 20
N = 200
t = linspace(0, tf, N)
U0 = array([1, 0, 0, 1])

U =  Cauchy(Kepler, t, U0, Euler) 

plt.plot(U[:,0] , U[:,1])
plt.show()

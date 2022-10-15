from numpy import array, linspace
import matplotlib.pyplot as plt
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv
from Functions.KeplerCode import Kepler
from ODEs.CauchyCode import Cauchy
from Richardson.ErrorCode import Error
from Richardson.ConvergenceCode import Conv

"Run"

#Defino tiempo final, numero divisiones y conodiciones iniciales
tf = 20
N = 100
t = linspace(0, tf, N)
U0 = array([1, 0, 0, 1])

#Calculo el error
E = Error(Kepler, t, U0, CN)
#Lo grafico
plt.plot(t, E[:,0])
plt.show()


#Calculo los logaritmos de los ratios de convergencia y numero de subdivisiones
[log_E, log_N] = Conv(Kepler, t, U0, CN, 8)
#Lo grafico
plt.plot(log_N, log_E)
plt.show()

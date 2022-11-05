from numpy import array, linspace, cos, sin
import matplotlib.pyplot as plt
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv, LF
from Functions.LinOscCode import LinOsc
from ODEs.CauchyCode import Cauchy
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv, LF
from Stability.RegionsCode import Regions
from Richardson.ErrorCode import Error
from Richardson.ConvergenceCode import Conv

#Defino tiempo final, numero divisiones y conodiciones iniciales
tf = 8
N = 80
t = linspace(0, tf, N)
U0 = array([1, 0])

#Integro el oscilador lineal
U = Cauchy(LinOsc, t, U0, RK4)
#Lo grafico
plt.plot(t, U[:,0],"r", label = "Solucion aproximada")
plt.plot(t, cos(t),"--",color = "r", label = "Solucion exacta")
plt.plot(t, U[:,1], "b", label = "Velocidad aproximada")
plt.plot(t, -sin(t),"--", color = "b", label = "Velocidad exacta")
plt.title("Oscilador lineal dt = 0.1, Leap-Frog")
plt.xlabel("t")
plt.ylabel("y(t)|dy(t)")
plt.legend(loc = "lower left")
plt.grid()
plt.show()


x = linspace(-5, 5, 100);
y = linspace(-5, 5, 100);

dt = array([0.01, 0.1, 0.5, 1])

stab_region = Regions(x,y,RK4)

CSF = plt.contourf(x, y, stab_region, levels = [0, 1],  colors=['#C0C0C0'] )
CS = plt.contour(x, y, stab_region, levels = [0.75, 1, 1.25] ) 
for t in range(len(dt)):
    plt.plot([0,0],[1*dt[t],-1*dt[t]], 'o', label = "Raíces del oscilador lineal con dt = " +str(dt[t]))
plt.clabel(CS, inline=1, fontsize=10)
plt.title('Regiones de estabilidad absoluta del Leap-Frog')
plt.xlabel("Re(|r|)")
plt.ylabel("Im(|r|)")
plt.legend(loc = "lower left")
plt.grid()
plt.show()








"""

#Calculo el error
E = Error(LinOsc, t, U0, LF)
#Lo grafico
plt.plot(t, E[:,0] ,"r", label = "Error de la posicion")
plt.plot(t, E[:,1] ,"b", label = "Error de la velocidad")
plt.title("Errores del Leap-Frog")
plt.xlabel("t")
plt.ylabel("Errores")
plt.legend(loc = "lower left")
plt.show()

#Calculo los logaritmos de los ratios de convergencia
[log_E, log_N] = Conv(LinOsc, t, U0, LF, 11)
#Lo grafico
plt.plot(log_N, log_E)
plt.title("Convergencia del Leap-Frog")
plt.xlabel("Log(N)")
plt.ylabel("Log(C)")
plt.show()

"""
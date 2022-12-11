from numpy import array, linspace, zeros, around
import matplotlib.pyplot as plt
from random import random

from ODEs.CauchyCode import Cauchy
from ODEs.SchemesCode import Euler, RK4, CN, EulerInv, LF
from ODEs.RKEmbCode import RKEmb
from Functions.Res3BodyCode import CR3BP, Lpoints, Stability

#Defino tiempo final, numero divisiones
N = int(1e4)
t = linspace(0, 1, N)
mu = 1.2151e-2 #Tierra-Luna
#mu = 3.0039e-7 #Sol-Tierra

# Creo una función que Cauchy pueda llamar ya que la nuestra tiene más inputs que (U, t)
def F(U,t):
   return CR3BP(U, mu)

#Calculo los puntos de Lagrange iniciando desde puntos cercanos
U0LP = array([[0.1, 0, 0, 0],[1.01, 0, 0, 0],[-0.1, 0, 0, 0],[0.8, 0.6, 0, 0],[0.8, -0.6, 0, 0]])
LagPoints = Lpoints(U0LP, 5, mu)
print(LagPoints)


#Genero una condiciones iniciales cercanas a un punto de lagrange
sel_LG = 5

U0 = zeros(4)
U0[0:2] = LagPoints[sel_LG-1,:]
ran = 1e-4*random()
U0 = U0 + ran

#Integro el problema circular restringido de los 3 cuerpos mediante un esquema temporal
U = Cauchy(F, t, U0, RKEmb)

#Estudio la estabilidad del punto
U0S = zeros(4)
U0S[0:2] = LagPoints[sel_LG-1,:]
eingvalues = Stability(U0S, mu)
print(around(eingvalues.real,8))

#Grafico
fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(U[:,0], U[:,1],'-',color = "r")
ax1.plot(-mu, 0, 'o', color = "g")
ax1.plot(1-mu, 0, 'o', color = "b")
for i in range(5):
    ax1.plot(LagPoints[i,0], LagPoints[i,1] , 'o', color = "k")
ax2.plot(U[:,0], U[:,1],'-',color = "r")
ax2.plot(LagPoints[sel_LG-1,0], LagPoints[sel_LG-1,1] , 'o', color = "k")
ax1.set_title("Orbital view")
ax2.set_title("Close-up")
fig.suptitle("Orbit around L2 with Embedded RK")
for ax in fig.get_axes():
    ax.set(xlabel='x', ylabel='y')
    ax.grid()

plt.show()

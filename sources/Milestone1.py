from numpy import array, zeros, transpose
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

"Nº pasos"
n = 200
"Tiempo final"
tf = 20
dt = tf/n
"C.I."
U0 = array([1,0,0,1]) 

"EULER"
U = zeros((4,n))
U[:,0] = U0
F = zeros(4)
for i in range(1, n):
    F = array([U[2,i-1],
               U[3,i-1],
               -U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2),
               -U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
    U[:,i] = U[:,i-1] + dt*F

plt.figure(1)
plt.plot(transpose(U[0,:]),transpose(U[1,:]))
plt.show()



"RK4"
U = zeros((4,n))
U[:,0] = U0
k1 = zeros(4)
k2 = zeros(4)
k3 = zeros(4)
k4 = zeros(4)
Unew1 = zeros(4)
Unew2 = zeros(4)
Unew3 = zeros(4)

for i in range(1, n):
    k1 = array([U[2,i-1],
                U[3,i-1],
                -U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2),
                -U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
    Unew1 = U[:,i-1] + dt*k1[:]/2
    k2 = array([Unew1[2],
                Unew1[3],
                -Unew1[0]/(Unew1[0]**2+Unew1[1]**2)**(3/2),
                -Unew1[1]/(Unew1[0]**2+Unew1[1]**2)**(3/2)])
    Unew2 = U[:,i-1] + dt*k2[:]/2
    k3 = array([Unew2[2],
                Unew2[3],
                -Unew2[0]/(Unew2[0]**2+Unew2[1]**2)**(3/2),
                -Unew2[1]/(Unew2[0]**2+Unew2[1]**2)**(3/2)])
    Unew3 = U[:,i-1] + dt*k3[:]
    k4 = array([Unew3[2],
                Unew3[3],
                -Unew3[0]/(Unew3[0]**2+Unew3[1]**2)**(3/2),
                -Unew3[1]/(Unew3[0]**2+Unew3[1]**2)**(3/2)])
    U[:,i] = U[:,i-1] + dt*(k1+2*k2+2*k3+k4)/6


plt.figure(2)
plt.plot(transpose(U[0,:]),transpose(U[1,:]))
plt.show()



"CRANK-NICHOLSON"

U = zeros((4,n))
U[:,0] = U0
F = zeros(4)

for i in range(1,n):
    F = array([U[2,i-1],
               U[3,i-1],
               -U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2),
               -U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
    def func_CN(x):
        return [x[0] - U[0,i-1] - (x[2] + F[0])*dt/2,
                x[1] - U[1,i-1] - (x[3] + F[1])*dt/2,
                x[2] - U[2,i-1] - (-x[0]/((x[0]**2+x[1]**2)**(3/2)) + F[2])*dt/2,
                x[3] - U[3,i-1] - (-x[1]/((x[0]**2+x[1]**2)**(3/2)) + F[3])*dt/2]
    U[:,i] = fsolve(func_CN, [U[:,i-1]])
        

plt.figure(3)
plt.plot(transpose(U[0,:]),transpose(U[1,:]))
plt.show()

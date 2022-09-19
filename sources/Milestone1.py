import numpy as np
import matplotlib.pyplot as plt

"No pasos"
n = 100
"Tiempo final"
tf = 10
dt = tf/n
"C.I."
U0 = np.array([1,0,0,1]) 




"EULER"
U = np.zeros((4,n))
U[:,0] = U0
F = np.zeros(4)
for i in range(1, n):
    F = np.array([U[2,i-1],U[3,i-1],-U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2),-U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
    U[:,i] = U[:,i-1] + dt*F

print(U)
plt.figure(1)
plt.plot(np.transpose(U[0,:]),np.transpose(U[1,:]))
plt.show()




"RK4"
U = np.zeros((4,n))
U[:,0] = U0
k1 = np.zeros(4)
k2 = np.zeros(4)
k3 = np.zeros(4)
k4 = np.zeros(4)
Unew1 = np.zeros(4)
Unew2 = np.zeros(4)
Unew3 = np.zeros(4)

for i in range(1, n):
    k1 = np.array([U[2,i-1],U[3,i-1],-U[0,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2),-U[1,i-1]/(U[0,i-1]**2+U[1,i-1]**2)**(3/2)])
    Unew1 = U[:,i-1] + dt*k1[:]/2
    k2 = np.array([Unew1[2],Unew1[3],-Unew1[0]/(Unew1[0]**2+Unew1[1]**2)**(3/2),-Unew1[1]/(Unew1[0]**2+Unew1[1]**2)**(3/2)])
    Unew2 = U[:,i-1] + dt*k2[:]/2
    k3 = np.array([Unew2[2],Unew2[3],-Unew2[0]/(Unew2[0]**2+Unew2[1]**2)**(3/2),-Unew2[1]/(Unew2[0]**2+Unew2[1]**2)**(3/2)])
    Unew3 = U[:,i-1] + dt*k3[:]
    k4 = np.array([Unew3[2],Unew3[3],-Unew3[0]/(Unew3[0]**2+Unew3[1]**2)**(3/2),-Unew3[1]/(Unew3[0]**2+Unew3[1]**2)**(3/2)])
    U[:,i] = U[:,i-1] + dt*(k1+2*k2+2*k3+k4)/6

print(U)
plt.figure(2)
plt.plot(np.transpose(U[0,:]),np.transpose(U[1,:]))
plt.show()




"Crank-Nicholson"











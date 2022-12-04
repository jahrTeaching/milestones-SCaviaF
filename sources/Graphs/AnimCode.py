from numpy import array, zeros
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

def Anim(Nb, r):
    def animate(num, data, line):
        line.set_data(data[0:2, :num])
        line.set_3d_properties(data[2, :num])
        return line
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    for j in range(Nb):
        x = zeros([N+1])
        y = zeros([N+1])
        z = zeros([N+1])
        c = colors[j]
        for i in range(N):
            x[i] = r[i,j,0]
            y[i] = r[i,j,1]
            z[i] = r[i,j,2]

        data = array([x, y, z])
        N = len(t)
        fig = plt.figure()
        ax = Axes3D(fig)

        line, = plt.plot(data[0], data[1], data[2], lw=1, c)
        line_ani = animation.FuncAnimation(fig, animate, frames=N, fargs=(data, line), interval=50, blit=False)

    plt.show()

    return

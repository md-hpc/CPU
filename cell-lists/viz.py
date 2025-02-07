from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
from math import floor
import os

DT, L = [float(x) for x in open("records/spec").read().split()]

N_PARTICLE = os.path.getsize("records/0") // 24

T = 1
while os.path.exists(f"records/{T}"):
    T += 1
fig = plt.figure()
ax = fig.add_subplot(projection='3d')

def gen(t, *fargs):
    particles = [[float(x) for x in line.split()] for line in open(f"records/{t}","r")]
    positions = np.zeros((3,N_PARTICLE))
    i = 0
    print(f"Drawing frame {t}")
    for i,p in enumerate(particles):
        positions[:,i] = np.array(p)
    return np.array_split(positions,3)    
        

def update(t, ax):
    ax.clear()
    xs, ys, zs = gen(t)
    ax.scatter(xs, ys, zs, c="r")
    ax.set_xlim3d([0.0, L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, L])
    ax.set_zlabel('Z')




# Setting the axes properties
ani = animation.FuncAnimation(fig, update, T, fargs=(ax,), interval=floor(DT*1000), repeat_delay=5000, blit=False)
ani.save('md.gif', writer='imagemagick')
plt.show()

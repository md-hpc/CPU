
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

import matplotlib.pyplot as plt
import numpy as np

timesteps = []

T = 0
L = 2.5 * 5
for line in open("particles"):
    line = [d.strip() for d in line.split(',')]
    t = int(line[0])
    l = max(float(d) for d in line[1:])
    T = max(t,T)

timesteps = [[] for _ in range(T+1)]
for line in open("particles"):
    line = [d.strip() for d in line.split(',')]
    t = int(line[0])
    r = np.array([float(d) for d in line[1:]])

    if np.any(r > L):
        print(f"dropping particle @ {t}: {r}")
    else:
        timesteps[t].append(r)

timesteps = [t for t in timesteps if len(t) > 0]
timesteps = [np.stack(t) for t in timesteps]

timesteps = timesteps[:200]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

def update(t, ax):
    print(f"Animating {t}")
    ax.clear()
    xs, ys, zs = np.split(timesteps[t],3,1)
    ax.scatter(xs, ys, zs, c="r", s=4)
    ax.set_xlim3d([0.0, L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0, L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0, L])
    ax.set_zlabel('Z')



# Setting the axes properties
ani = animation.FuncAnimation(fig, update, len(timesteps), fargs=(ax,), interval=10, repeat_delay=5000, blit=False)
ani.save('md.gif', writer='imagemagick')
plt.show()

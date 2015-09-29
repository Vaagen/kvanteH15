import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import time
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

t = time.clock()

def setV(V):
    return V

def setPsi(Nx, Ny, Nt):
    psi_r = np.zeros((Nx, Ny, Nt))
    psi_r
    return psi_r, psi_i

# set constants
Lx = 1e-3
Ly = 1e-3
Nx = 1000
Ny= 1000
dx = Lx / Nx
dy = Ly / Ny
dt = min([dx, dy])**2 * 1/4
t = 0
Nt = 100
hbar = 1.0545718e-34
m = 10e-30
V = np.zeros((Nx, Ny))
V = setV(V)

psi_r = np.zeros((Nx, Ny, Nt))
psi_i = np.zeros((Nx, Ny, Nt))
P = np.zeros(Nt)

c1x = hbar * dt / 2 / m / dx**2
c1y = hbar * dt / 2 / m / dy**2
c2V = V * dt / hbar + 2 * c1x + 2 * c1y
i = 0

fig = plt.figure()
ax = fig.gca(projection='3d')
X = np.arange(1, Nx-1, 1)
Y = np.arange(1, Ny-1, 1)
X, Y = np.meshgrid(X, Y)
while i < Nt - 1:
    psi_r[1:Nx-1, 1:Ny-1, i + 1] = psi_r[1:Nx-1, 1:Ny-1, i] - \
            c1x * (psi_i[2:Nx, 1:Ny-1, i] + psi_i[0:Nx-2, 1:Ny-1, i]) - \
            c1y * (psi_i[1:Nx-1, 2:Ny, i] + psi_i[1:Nx-1, 0:Ny-2, i]) + \
            c2V[1:Nx-1, 1:Ny-1] * psi_i[1:Nx-1, 1:Ny-1, i]
    psi_i[1:Nx-1, 1:Ny-1, i + 1] = psi_i[1:Nx-1, 1:Ny-1, i] + \
            c1x * (psi_r[2:Nx, 1:Ny-1, i] + psi_r[0:Nx-2, 1:Ny-1, i]) + \
            c1y * (psi_r[1:Nx-1, 2:Ny, i] + psi_r[1:Nx-1, 0:Ny-2, i]) - \
            c2V[1:Nx-1, 1:Ny-1] * psi_r[1:Nx-1, 1:Ny-1, i]
    if i % 50 == 0:
        surf = ax.plot_surface(X, Y, psi_r[1:Nx-1, 1:Ny-1, i + 1], rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        ax.set_zlim(-1.01, 1.01)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

        fig.colorbar(surf, shrink=0.5, aspect=5)

        plt.draw()
    print i
    i += 1
T = time.clock() - t
min = round(T / 60)
sec = T % 60
print "min " + str(min) + " sec " + str(sec)
plt.show()


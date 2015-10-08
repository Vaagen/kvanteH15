# this script plot simulations run by 'Schrodinger.h' and 'Schrodinger.cpp'
# one only needs to specify the general filename common for all files produced in 'Schrodinger.h' and 'Schrodinger.cpp' (the filename argument in the fun member function)
# as well as placing these files in the same directory as this file, 'plotSchroinger.py', or in a subdirectory of this directory

fileName = "test_free_electron"

# one should not be needing to do changes to the rest of the script

import numpy as np
import matplotlib
matplotlib.use('TKAgg') # this is done to allow blit = True in FuncAnimation on mac
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import os
import time

# find location of file, 'name'
def find(name):
    for root, dirs, files in os.walk(os.path.dirname(os.path.realpath(__file__))):
        if name in files:
            return os.path.join(root, name)

# load variables used in simulation
variableFile = open(find(fileName + "_variables.txt"), mode = 'r')
line = variableFile.readline()
numOfDim = int(line)
line = variableFile.readline()
Lx1 = float(line)
line = variableFile.readline()
Lx2 = float(line)
line = variableFile.readline()
Lx3 = float(line)
line = variableFile.readline()
Nx1 = int(line)
line = variableFile.readline()
Nx2 = int(line)
line = variableFile.readline()
Nx3 = int(line)
line = variableFile.readline()
Nt = int(line)
line = variableFile.readline()
dx1 = float(line)
line = variableFile.readline()
dx2 = float(line)
line = variableFile.readline()
dx3 = float(line)
line = variableFile.readline()
dt = float(line)
line = variableFile.readline()
m = float(line)
line = variableFile.readline()
p = float(line)
line = variableFile.readline()
k = float(line)
line = variableFile.readline()
startX1 = int(line)
line = variableFile.readline()
startX2 = int(line)
line = variableFile.readline()
startX3 = int(line)
line = variableFile.readline()
V0 = float(line)
line = variableFile.readline()
VThickness = float(line)
line = variableFile.readline()
situation = int(line)
line = variableFile.readline()
potential = int(line)
line = variableFile.readline()
probDistrb = int(line)
line = variableFile.readline()
SDx1 = float(line)
line = variableFile.readline()
SDx2 = float(line)
line = variableFile.readline()
SDx3 = float(line)
line = variableFile.readline()
plotDensityX1 = int(line)
line = variableFile.readline()
plotDensityX2 = int(line)
line = variableFile.readline()
plotDensityX3 = int(line)
line = variableFile.readline()
plotDensityT = int(line)
variableFile.close()

# get data from plotFile
dt = np.dtype("f8")
plotProbabilityFile = np.fromfile(find(fileName + "_plot_probability"), dtype=dt)
plotPsiRFile = np.fromfile(find(fileName + "_plot_psi_r"), dtype=dt)
plotPsiIFile = np.fromfile(find(fileName + "_plot_psi_i"), dtype=dt)

'''
x = np.linspace(0,Lx1,Nx1/plotDensityX1)
plt.plot(x, plotPsiRFile[0:Nx1/plotDensityX1], 'r.')
plt.plot(x, plotPsiRFile[100*Nx1/plotDensityX1:(100+1)*Nx1/plotDensityX1], 'g.')
plt.show()
'''

# for i in range(0,Nt/plotDensityT):
#    print plotProbabilityFile[i*Nx1/plotDensityX1 + 100]

fig = plt.figure()
ax = plt.axes(xlim=(0, Lx1), ylim=(-1, 1))
probPlot, = ax.plot([], [], 'k', lw = 1, label = 'Probability')
psiRPlot, = ax.plot([], [], 'b', lw = 1, label = 'Real part') # only used for 1D
psiIPlot, = ax.plot([], [], 'r', lw = 1, label = 'Imaginary part') # only used for 1D
plt.legend(loc = 'lower right')
if numOfDim == 2:
    probPlot, = ax.contourf([], [], [])

# initialization function: plot the background of each frame
def init1D():
    probPlot.set_data([], [],)
    return probPlot

def init2D():
    #plotData = ax.contourf([], [], [], 500)
    return

def init3D():
    return

x1 = np.linspace(0,Lx1,Nx1/plotDensityX1)
x2 = np.linspace(0,Lx2,Nx2/plotDensityX2)
x3 = np.linspace(0,Lx3,Nx3/plotDensityX3)
# animation function.  This is called sequentially
def animate1D(i):
    probPlot.set_data(x1, plotProbabilityFile[Nx1/plotDensityX1*i:Nx1/plotDensityX1*(i+1)])
    psiRPlot.set_data(x1, plotPsiRFile[Nx1/plotDensityX1*i:Nx1/plotDensityX1*(i+1)])
    psiIPlot.set_data(x1, plotPsiIFile[Nx1/plotDensityX1*i:Nx1/plotDensityX1*(i+1)])
    return probPlot, psiRPlot, psiIPlot

def animate2D(i):
    
    return

def animate3D(i):
    return


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = None
if numOfDim == 1:
    anim = animation.FuncAnimation(fig, animate1D, frames = Nt/plotDensityT, interval=20, blit=True, repeat = False)
elif numOfDim == 2:
    anim = animation.FuncAnimation(fig, animate2D, init_func=init2D, frames = Nt/plotDensityT, interval=20, blit=True, repeat = False)
elif numOfDim == 3:
    anim = animation.FuncAnimation(fig, animate3D, init_func=init3D, frames = Nt/plotDensityT, interval=20, blit=True, repeat = False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show(block = False)
plt.close()

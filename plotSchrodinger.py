# this script plot simulations run by 'Schrodinger.h' and 'Schrodinger.cpp'
# one only needs to specify the general filename common for all files produced in 'Schrodinger.h' and 'Schrodinger.cpp' (the filename argument in the fun member function)

fileName = "test_free_electron"

# one should not be needing to do changes to the rest of the script

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import animation
import os

# find location of file, 'name'
def find(name):
    for root, dirs, files in os.walk(os.path.dirname(os.path.realpath(__file__))):
        if name in files:
            return os.path.join(root, name)

# load variables used in simulation
variableFile = open(find(fileName + "_variables.txt"), mode = 'r')
line = variableFile.readline()
numOfDim = float(line)
line = variableFile.readline()
Lx1 = float(line)
line = variableFile.readline()
Lx2 = float(line)
line = variableFile.readline()
Lx3 = float(line)
line = variableFile.readline()
Nx1 = float(line)
line = variableFile.readline()
Nx2 = float(line)
line = variableFile.readline()
Nx3 = float(line)
line = variableFile.readline()
Nt = float(line)
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
startX1 = float(line)
line = variableFile.readline()
startX2 = float(line)
line = variableFile.readline()
startX3 = float(line)
line = variableFile.readline()
V0 = float(line)
line = variableFile.readline()
VThickness = float(line)
line = variableFile.readline()
situation = float(line)
line = variableFile.readline()
potential = float(line)
line = variableFile.readline()
probDistrb = float(line)
line = variableFile.readline()
SDx1 = float(line)
line = variableFile.readline()
SDx2 = float(line)
line = variableFile.readline()
SDx3 = float(line)
line = variableFile.readline()
plotDensityX1 = float(line)
line = variableFile.readline()
plotDensityX2 = float(line)
line = variableFile.readline()
plotDensityX3 = float(line)
line = variableFile.readline()
plotDensityT = float(line)
variableFile.close()

# get data from plotFile
dt = np.dtype("f8")
plotFile = np.fromfile(find(fileName + "_plot"), dtype=dt)
if numOfDim == 1:
    plotFile = np.split(plotFile, 3) # only simulations for 1 dim store psi_r and psi_i in plotFile

fig = plt.figure()
ax = plt.subplot(1,1,1)
plotData
if numOfDim == 1:
    plotData = ax.plot([],[],[],[],[],[])
else if numOfDim == 2:
    plotData = ax.contourf([], [], [], 500)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, Lx1, Nx1)
    y = plotFile[0][i]
    line.set_data(x, y)
    return line,



# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = fig.add_subplot(111)
line, = ax.plot([], [])

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init, interval=20, blit=False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

plt.show()

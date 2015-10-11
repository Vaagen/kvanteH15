# this script plot simulations run by 'Schrodinger.h' and 'Schrodinger.cpp'
# one only needs to specify the general filename common for all files produced in 'Schrodinger.h' and 'Schrodinger.cpp' (the filename argument in the fun member function)
# as well as placing these files in the same directory as this file, 'plotSchroinger.py', or in a subdirectory of this directory

fileName = "test_free_electron"
animationType2D = "frameByFrame"  #choose between "animation" or "frameByFrame"

# one should not be needing to do changes to the rest of the script

import numpy as np
import matplotlib
matplotlib.use('TKAgg') # this is done to allow blit = True in FuncAnimation on mac
# import proper graphics back-end for Mac OS X
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import pylab
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
Vmax = float(line)
line = variableFile.readline()
startEnergy = float(line)
line = variableFile.readline()
finalEnergy = float(line)
line = variableFile.readline()
finalProb = float(line)
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
plotSpacingX1 = int(line)
line = variableFile.readline()
plotSpacingX2 = int(line)
line = variableFile.readline()
plotSpacingX3 = int(line)
line = variableFile.readline()
plotSpacingT = int(line)
variableFile.close()

print "The final energy is ",
print finalEnergy / startEnergy,
print " times the start energy."
print "The final probability of finding the particle is: ",
print finalProb

# get data from plotFile
dt = np.dtype("f8")
plotProbabilityFile = np.fromfile(find(fileName + "_plot_probability"), dtype=dt)
plotPsiRFile = np.fromfile(find(fileName + "_plot_psi_r"), dtype=dt)
plotPsiIFile = np.fromfile(find(fileName + "_plot_psi_i"), dtype=dt)
potentialFile = np.fromfile(find(fileName + "_potential"), dtype=dt)


'''
x = np.linspace(0,Lx1,Nx1/plotSpacingX1)
plt.plot(x, plotPsiRFile[0:Nx1/plotSpacingX1], 'r.')
plt.plot(x, plotPsiRFile[100*Nx1/plotSpacingX1:(100+1)*Nx1/plotSpacingX1], 'g.')
plt.show()
'''

startTime = time.clock()

fig = plt.figure()
'''
scaleConstPsi = 0.8 * max(plotProbabilityFile) / (max([max(plotPsiRFile),max(plotPsiIFile)]))
ax = plt.axes(xlim=(0, Lx1), ylim=(1.1 * scaleConstPsi * min([min(plotPsiRFile), min(plotPsiIFile)]), 1.1 * max(plotProbabilityFile)))
probPlot, = ax.plot([], [], 'k', lw = 1, label = 'Probability')
psiRPlot, = ax.plot([], [], 'b', lw = 1, label = 'Real part') # only used for 1D
psiIPlot, = ax.plot([], [], 'r', lw = 1, label = 'Imaginary part') # only used for 1D
plt.legend(loc = 'lower right')
'''

x1 = np.linspace(0,Lx1,Nx1/plotSpacingX1)
x2 = np.linspace(0,Lx2,Nx2/plotSpacingX2)
x3 = np.linspace(0,Lx3,Nx3/plotSpacingX3)

scaleConstEnergy = 0.5 * max(plotProbabilityFile) / startEnergy

ax = fig.gca(projection = '3d')
if numOfDim == 1:
    plt.plot(x1, scaleConstEnergy * potentialFile, ':k', zorder=0)
    pylab.fill(x1, scaleConstEnergy * potentialFile, facecolor='y', alpha=0.2, zorder=0)
'''
if numOfDim == 2:
    probPlot = ax3d.plot_surface([], [], [], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    x1, x2 = np.meshgrid(x1, x2)
'''
z = plotProbabilityFile[0:Nx1/plotSpacingX1*Nx2/plotSpacingX2].reshape(Nx2/plotSpacingX2,Nx1/plotSpacingX1)
x1, x2 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx2,Nx2/plotSpacingX2))
probPlot = ax.plot_surface(x1, x2, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False);

print "There are ",
print Nt/plotSpacingT,
print " number of frames."
print "press enter to go to next frame"
if numOfDim == 2 and animationType2D == "frameByFrame":
    for t in range(0,Nt/plotSpacingT):
        plt.cla()
        z = plotProbabilityFile[Nx1/plotSpacingX1*Nx2/plotSpacingX2*t:Nx1/plotSpacingX1*Nx2/plotSpacingX2*(t+1)].reshape(Nx2/plotSpacingX2,Nx1/plotSpacingX1)
        probPlot = ax.plot_surface(x1, x2, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
        plt.show(block = False)
        raw_input()

# initialization function: plot the background of each frame
def init1D():
    energy = np.linspace(scaleConstEnergy * startEnergy, scaleConstEnergy * startEnergy, Nx1/plotSpacingX1)
    energyPlot, = ax.plot(x1, energy)
    probPlot.set_data([], [])
    psiRPlot.set_data([], [])
    psiIPlot.set_data([], [])
    return probPlot, psiRPlot, psiIPlot, energyPlot,

def init2D():
    probPlot = ax.plot_surface([], [], [], rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx2,Nx2/plotSpacingX2))
    return probPlot,

def init3D():
    return


# animation function.  This is called sequentially
def animate1D(i):
    probPlot.set_data(x1, plotProbabilityFile[Nx1/plotSpacingX1*i:Nx1/plotSpacingX1*(i+1)])
    psiRPlot.set_data(x1, scaleConstPsi * plotPsiRFile[Nx1/plotSpacingX1*i:Nx1/plotSpacingX1*(i+1)])
    psiIPlot.set_data(x1, scaleConstPsi * plotPsiIFile[Nx1/plotSpacingX1*i:Nx1/plotSpacingX1*(i+1)])
    return probPlot, psiRPlot, psiIPlot,

def animate2D(i, x1, x2, probPlot):
    z = plotProbabilityFile[Nx1/plotSpacingX1*Nx2/plotSpacingX2*i:Nx1/plotSpacingX1*Nx2/plotSpacingX2*(i+1)].reshape(Nx2/plotSpacingX2,Nx1/plotSpacingX1)
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx2,Nx2/plotSpacingX2))
    probPlot = ax.plot_surface(x1, x2, z, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    return probPlot,

def animate3D(i):
    return


# call the animator.  blit=True means only re-draw the parts that have changed.
anim = None
print Nt/plotSpacingT * 50 / 1000
if numOfDim == 1:
    anim = animation.FuncAnimation(fig, animate1D, init_func=init1D, frames = Nt/plotSpacingT, interval=20, blit=True, repeat = False)
elif numOfDim == 2 and animationType2D == "animation":
    anim = animation.FuncAnimation(fig, animate2D, frames = Nt/plotSpacingT, interval=50, blit=False, repeat = False, fargs=(x1,x2,probPlot))
elif numOfDim == 3:
    anim = animation.FuncAnimation(fig, animate3D, init_func=init3D, frames = Nt/plotSpacingT, interval=20, blit=True, repeat = False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
#anim.save('basic_animation.mp4', fps=30, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])
if numOfDim != 2 or animationType2D == "animation":
    plt.show(block = False) # the 'block = False' somehow makes the animation unstable and make it randomly quit prematurly

plt.close()
print "Seconds used to run animation: ",
print time.clock() - startTime
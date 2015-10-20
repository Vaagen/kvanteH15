# this script plot simulations run by 'Schrodinger.h' and 'Schrodinger.cpp'
# one only needs to specify the general filename common for all files produced in 'Schrodinger.h' and 'Schrodinger.cpp' (the filename argument in the fun member function)
# as well as placing these files in the same directory as this file, 'plotSchroinger.py', or in a subdirectory of this directory

fileName = "test_free_electron"

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
import math
import customColormaps as ccmaps

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

# get data from xFile'ene
dt = np.dtype("f8")
plotProbabilityFile = np.fromfile(find(fileName + "_plot_probability"), dtype=dt)
plotPsiRFile = np.fromfile(find(fileName + "_plot_psi_r"), dtype=dt)
plotPsiIFile = np.fromfile(find(fileName + "_plot_psi_i"), dtype=dt)
potentialFile = np.fromfile(find(fileName + "_potential"), dtype=dt)

fig = plt.figure(1)

# axes used for simulation in 3D
ax1 = None
ax2 = None
ax3 = None
ax4 = None

scaleConstPsi = 0.8 * max(plotProbabilityFile) / (max([max(plotPsiRFile),max(plotPsiIFile)]))
ax = plt.axes(xlim=(0, Lx1), ylim=(1.1 * scaleConstPsi * min([min(plotPsiRFile), min(plotPsiIFile)]), 1.1 * max(plotProbabilityFile)))
probPlot, = ax.plot([], [], 'k', lw = 1, zorder = 3, label = 'Probability')
psiRPlot, = ax.plot([], [], 'b', lw = 1, zorder = 2, label = 'Re($\Psi$)') # only used for 1D
psiIPlot, = ax.plot([], [], 'r', lw = 1, zorder = 2, label = 'Im($\Psi$)') # only used for 1D


x1 = np.linspace(0,Lx1,Nx1/plotSpacingX1)
x2 = np.linspace(0,Lx2,Nx2/plotSpacingX2)
x3 = np.linspace(0,Lx3,Nx3/plotSpacingX3)

scaleConstEnergy = 0.5 * max(plotProbabilityFile) / startEnergy

# initialization function: plot the background of each frame
def init1D():
    plt.plot(x1, scaleConstEnergy * potentialFile, ':k', zorder=0)
    pylab.fill(x1, scaleConstEnergy * potentialFile, facecolor='y', alpha=0.2, zorder=0, label = 'Potential')
    energy = np.linspace(scaleConstEnergy * startEnergy, scaleConstEnergy * startEnergy, Nx1/plotSpacingX1)
    energyPlot, = ax.plot(x1, energy, 'g', zorder = 1, label = "Energy")
    plt.legend(loc = 'lower right', prop={'size':10})
    probPlot.set_data([], [])
    psiRPlot.set_data([], [])
    psiIPlot.set_data([], [])
    return probPlot, psiRPlot, psiIPlot, energyPlot,

def init2D():
    plt.cla()
    z = potentialFile[0:Nx1*Nx2:plotSpacingX1].reshape(Nx2,Nx1/plotSpacingX1)[0:Nx2:plotSpacingX2,:]
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx2,Nx2/plotSpacingX2))
    colorLimit = max(abs(potentialFile))
    if colorLimit == 0:
        colorLimit = 1
    potentialPlot = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('white_black'), zorder = 2, vmin = - colorLimit, vmax = colorLimit)
    return potentialPlot

def init3D():
    ax1 = plt.add_subplot(221)
    ax2 = plt.add_subplot(222)
    ax3 = plt.add_subplot(223)
    ax4 = plt.add_subplot(224)
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
    probPlot = plt.contourf(x1, x2, z, cmap=ccmaps.cmap('vaagen_colorscale'), zorder = 1, vmin = 0, vmax = max(plotProbabilityFile))
    return probPlot,

def animate3D(i):
    probMat = plotProbabilityFile[Nx1/plotSpacingX1*Nx2/plotSpacingX2*Nx3/plotSpacingX3*i:Nx1/plotSpacingX1*Nx2/plotSpacingX2*Nx3/plotSpacingX3*(i+1)].reshape(Nx3/plotSpacingX3,Nx2/plotSpacingX2,Nx1/plotSpacingX1)
    z11 = np.sum(probMat, axis=0)
    x1, x2 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx2,Nx2/plotSpacingX2))
    plt.sca(ax1)
    probPlot11 = plt.contourf(x1, x2, z11, cmap=plt.cm.gnuplot2_r, alpha = 1, zorder = 1, label = 'xy-axis')
    z12 = np.sum(probMat, axis=1)
    x1, x3 = np.meshgrid(np.linspace(0,Lx1,Nx1/plotSpacingX1), np.linspace(0,Lx3,Nx3/plotSpacingX3))
    plt.sca(ax2)
    probPlot12 = plt.contourf(x1, x3, z12, cmap=plt.cm.gnuplot2_r, alpha = 1, zorder = 1, label = 'xz-axis')
    z21 = np.sum(probMat, axis=2)
    x2, x3 = np.meshgrid(np.linspace(0,Lx2,Nx2/plotSpacingX2), np.linspace(0,Lx3,Nx3/plotSpacingX3),)
    plt.sca(ax3)
    probPlot21 = plt.contourf(x2, x3, z21, cmap=plt.cm.autumn_r, alpha = 1, zorder = 1, label = 'yz-axis')
    # http://stackoverflow.com/questions/9164950/python-matplotlib-plot3d-with-a-color-for-4d
    return


# call the animator.  blit=True means only re-draw the parts that have changed.
startTime = time.clock()
anim = None
animationTime = 10; #sec
if numOfDim == 1:
    anim = animation.FuncAnimation(fig, animate1D, init_func=init1D, frames = Nt/plotSpacingT, interval=20, blit=True, repeat = False)
elif numOfDim == 2:
    anim = animation.FuncAnimation(fig, animate2D, init_func=init2D, frames = Nt/plotSpacingT, interval=1000*animationTime/Nt*plotSpacingT, blit=False, repeat = False, fargs=(x1,x2,probPlot))
elif numOfDim == 3:
    anim = animation.FuncAnimation(fig, animate3D, init_func=init3D, frames = Nt/plotSpacingT, interval=20, blit=True, repeat = False)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html
anim.save('animMovies/basic_animation.mp4', fps=Nt/plotSpacingT/animationTime, extra_args=['-vcodec', 'libx264', '-pix_fmt', 'yuv420p'])

print "Seconds used to run animation: ",
print time.clock() - startTime
//
//  Schrodinger2D.cpp
//  schrodingerFD2D
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger2D.h"



double Lx = 1 * pow(10, -3);
double Ly = 1 * pow(10, -3);
int Nx = 1000;
int Ny = 1000;
double dx = Lx / Nx;
double dy = Ly / Ny;
double dt = dx * dx * 1/4;
if (dy < dx){
    dt = dy * dy * 1/4;
}
int Nt = 1000;
double hbar = 1.0545718 * pow(10, -34);
}
/*
 Lx = 1e-3
 Ly = 1e-3
 Nx = 1000
 Ny= 1000
 dx = Lx / Nx
 dy = Ly / Ny
 dt = min([dx, dy])**2 * 1/4
 t = 0
 Nt = 1000
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
 
 nthplt = 10
 X = np.arange(1, Nx-1, nthplt)
 Y = np.arange(1, Ny-1, nthplt)
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
 if i % 50 == 49:
 surf = ax.plot_surface(X, Y, psi_r[1:Nx-1:nthplt, 1:Ny-1:nthplt, i + 1], rstride=1, cstride=1, cmap=cm.coolwarm,
 linewidth=0, antialiased=False)
 plt.draw()
 print i
 i += 1
 */
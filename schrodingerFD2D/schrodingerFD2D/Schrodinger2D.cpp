//
//  Schrodinger2D.cpp
//  schrodingerFD2D
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger2D.h"
#include <string>

using namespace std;
//PUBLIC MEMBER FUNCTIONS
void Schrodinger2D::run(){
    
}

void Schrodinger2D::makeInitCondition(){
    
}

void Schrodinger2D::contSim(string filename, unsigned int numOfTimesteps){
    
}

//PRIVATE MEMBER FUNCTIONS
void Schrodinger2D::setV(Potential type){
    switch (type) {
        case FREE:
            for (int x = 0; x < Nx; x++){
                for (int y = 0; y < Ny; y++){
                    V[Nx*y+x] = 0;
                }
            }
            break;
        default:
            break;
    }
}

void Schrodinger2D::initState(string type){
    
}

void Schrodinger2D::finiteDifference(){
    
}
/*
 
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
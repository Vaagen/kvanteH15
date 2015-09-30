//
//  Schrodinger2D.cpp
//  schrodingerFD2D
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger2D.h"

using namespace std;
//PUBLIC MEMBER FUNCTIONS
void Schrodinger2D::run(){
    int choice = -1;
    while (choice < 0 || choice >= situationString.size()){
        cout << "Choose a number correspolding to the situation you would like to simulate:" << endl;
        for (int i = 0; i < situationString.size(); i++){
            cout << i + 1 << ". " << situationString[i] << endl;
        }
        cin >> choice;
    }
    runSituation(static_cast<Situation>(choice));
    
}

void Schrodinger2D::runSituation(Situation situation){
    // standard settings (you should generly override these in you 'Situation' unless it is benifitial to change them in all situations)
    Lx = 0.001;
    Ly = 0.001;
    Nx = 1000;
    Ny = 1000;
    Nt = 1000;
    dx = Lx / Nx;
    dy = Ly / Ny;
    dt = 1/4 * pow((Lx / Nx),2);
    if (Ly / Ny < Lx / Nx){
        dt = 1/4 * pow((Ly / Ny),2);
    }
    startX = Lx / 4;
    startY = Ly / 2;
    plotDensityX = 1;
    plotDensityY = 1;
    plotDensityT = 50;
    
    switch (situation) {
        case FREE_ELECTRON:
            m = pow(10, -30);
            potential = FREE;
            probDistrb = GAUSSIAN;
            SDx = 1;
            SDy = 1;
            p = 10;
            break;
        default:
            break;
    }
    setV();
    makeInitState();
    finiteDifference();
}

void Schrodinger2D::contSim(string filename, unsigned int numOfTimesteps){
    
}

//PRIVATE MEMBER FUNCTIONS
void Schrodinger2D::setV(){
    switch (potential) {
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

void Schrodinger2D::makeInitState(){
    psi_r = new double [Nx * Ny * Nt];
    psi_i = new double [Nx * Ny * Nt];
    switch (probDistrb) {
        case GAUSSIAN:
            for (int i = 0; i < Nx; i++){
                for (int j = 0; j < Ny; j++){
                    psi_r[i * Ny + j] = exp(-pow(dx * i - startX, 2) / (2 * SDx) - pow(dy * j - startY, 2) / (2 * SDy)) * cos(p * (dx * i));
                    psi_i[i * Ny + j] = exp(-pow(dx * i - startX, 2) / (2 * SDx) - pow(dy * j - startY, 2) / (2 * SDy)) * sin(p * (dx * i));
                }
            }
            break;
            
        default:
            break;
    }
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
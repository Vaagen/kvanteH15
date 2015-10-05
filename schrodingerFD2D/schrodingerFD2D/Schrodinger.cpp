//
//  Schrodinger.cpp
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger.h"

using namespace std;

//PUBLIC MEMBER FUNCTIONS
void Schrodinger::run(Situation situation, string filename){
    // standard settings (you should generely override these in you 'Situation', unless it is benifitial to change them in all situations)
    numOfDim = 1;
    Lx1 = 0.001;
    Lx2 = 0.001;
    Lx3 = 0.001;
    Nx1 = 1000;
    Nx2 = 1000;
    Nx3 = 1000;
    if (numOfDim == 1){
        Nx2 = 1;
        Nx3 = 1;
    } else if (numOfDim == 2){
        Nx3 = 1;
    }
    Nt = 1000;
    dx1 = Lx1 / Nx1;
    dx2 = Lx2 / Nx2;
    dx3 = Lx3 / Nx3;
    dt = 1/4 * pow((Lx1 / Nx1),2);
    if (Lx2 / Nx2 < Lx1 / Nx1){
        dt = 1/4 * pow((Lx2 / Nx2),2);
    }
    startX1 = Lx1 / 4;
    startX2 = Lx2 / 2;
    plotDensityX1 = 1;
    plotDensityX2 = 1;
    plotDensityT = 50;
    
    switch (situation) {
        case FREE_ELECTRON_1D:
            numOfDim = 1;
            m = pow(10, -30);
            potential = FREE;
            probDistrb = GAUSSIAN_1D;
            SDx1 = 1;
            SDx2 = 1;
            p = 10;
            break;
        default:
            break;
    }
    setV();
    makeInitState();
    finiteDifference();
}

void Schrodinger::contSim(string filename, unsigned int numOfTimesteps){
}

Schrodinger::Schrodinger(){
    V = nullptr;
    psi_r = nullptr;
    psi_i = nullptr;
}

Schrodinger::~Schrodinger(){
    delete [] V;
    delete [] psi_r;
    delete [] psi_i;
}

//PRIVATE MEMBER FUNCTIONS
void Schrodinger::setV(){
    //ofstream ofs( "atest.txt", ios::binary );
    V = new double [Nx1 * Nx2 * Nx3];
    switch (potential) {
        case FREE:
            for (int x1 = 0; x1 < Nx1; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                    for (int x3 = 0; x3 < Nx3; x3++){
                        V[Nx1*Nx2*x3+Nx1*x2+x1] = 0;
                    }
                }
            }
            break;
        default:
            break;
    }
}

void Schrodinger::makeInitState(){
    psi_r = new double [Nx1 * Nx2 * Nt];
    psi_i = new double [Nx1 * Nx2 * Nt];
    switch (probDistrb) {
        case GAUSSIAN_2D:
            for (int i = 0; i < Nx1; i++){
                for (int j = 0; j < Nx2; j++){
                    psi_r[i * Nx2 + j] = exp(-pow(dx1 * i - startX1, 2) / (2 * SDx1) - pow(dx1 * j - startX2, 2) / (2 * SDx2)) * cos(p * (dx1 * i));
                    psi_i[i * Nx2 + j] = exp(-pow(dx1 * i - startX1, 2) / (2 * SDx1) - pow(dx2 * j - startX2, 2) / (2 * SDx2)) * sin(p * (dx1 * i));
                }
            }
            break;
            
        default:
            break;
    }
}

void Schrodinger::finiteDifference(){
    if (numOfDim == 1){
        finiteDifference1D();
    } else if (numOfDim == 2){
        finiteDifference2D();
    } else {
        finiteDifference3D();
    }
}

void Schrodinger::finiteDifference1D(){
    
}

void Schrodinger::finiteDifference2D(){
    
}

void Schrodinger::finiteDifference3D(){
    
}

/*
 
 psi_r = np.zeros((Nx1, Nx12, Nt))
 psi_i = np.zeros((Nx1, Nx12, Nt))
 P = np.zeros(Nt)
 
 c1x = hbar * dt / 2 / m / dx**2
 c1y = hbar * dt / 2 / m / dy**2
 c2V = V * dt / hbar + 2 * c1x + 2 * c1y
 i = 0
 
 fig = plt.figure()
 ax = fig.gca(projection='3d')
 
 nthplt = 10
 X = np.arange(1, Nx1-1, nthplt)
 Y = np.arange(1, Nx12-1, nthplt)
 X, Y = np.meshgrid(X, Y)
 while i < Nt - 1:
 psi_r[1:Nx1-1, 1:Nx12-1, i + 1] = psi_r[1:Nx1-1, 1:Nx12-1, i] - \
 c1x * (psi_i[2:Nx1, 1:Nx12-1, i] + psi_i[0:Nx1-2, 1:Nx12-1, i]) - \
 c1y * (psi_i[1:Nx1-1, 2:Nx12, i] + psi_i[1:Nx1-1, 0:Nx12-2, i]) + \
 c2V[1:Nx1-1, 1:Nx12-1] * psi_i[1:Nx1-1, 1:Nx12-1, i]
 psi_i[1:Nx1-1, 1:Nx12-1, i + 1] = psi_i[1:Nx1-1, 1:Nx12-1, i] + \
 c1x * (psi_r[2:Nx1, 1:Nx12-1, i] + psi_r[0:Nx1-2, 1:Nx12-1, i]) + \
 c1y * (psi_r[1:Nx1-1, 2:Nx12, i] + psi_r[1:Nx1-1, 0:Nx12-2, i]) - \
 c2V[1:Nx1-1, 1:Nx12-1] * psi_r[1:Nx1-1, 1:Nx12-1, i]
 if i % 50 == 49:
 surf = ax.plot_surface(X, Y, psi_r[1:Nx1-1:nthplt, 1:Nx12-1:nthplt, i + 1], rstride=1, cstride=1, cmap=cm.coolwarm,
 linewidth=0, antialiased=False)
 plt.draw()
 print i
 i += 1
 */
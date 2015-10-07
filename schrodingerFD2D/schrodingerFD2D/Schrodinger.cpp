//
//  Schrodinger.cpp
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include "Schrodinger.h"
#include <fstream>
#include "float.h"

using namespace std;

//PUBLIC MEMBER FUNCTIONS
void Schrodinger::run(Situation situation, string filename){
    // standard settings (you should generely override these in you 'Situation', unless it is benifitial to change them in all situations)
    this->filename = filename;
    numOfDim = 1;
    Lx1 = 0.001;
    Lx2 = 0.001;
    Lx3 = 0.001;
    Nx1 = 1000;
    Nx2 = 1000;
    Nx3 = 1000;
    Nt = 1000;
    V0 = 1;
    VThickness = 1;
    
    SDx1 = 1;
    SDx2 = 1;
    SDx3 = 1;
    
    plotDensityX1 = 1;
    plotDensityX2 = 1;
    plotDensityX3 = 1;
    plotDensityT = 50;
    
    // set the situation
    // this is where you add new situations (along with adding i new situation in the enum Situation)
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
    
    if (numOfDim == 1){
        Nx2 = 1;
        Nx3 = 1;
    } else if (numOfDim == 2){
        Nx3 = 1;
    }
    dx1 = Lx1 / Nx1;
    dx2 = Lx2 / Nx2;
    dx3 = Lx3 / Nx3;
    dt = 1/4 * pow((Lx1 / Nx1),2); // should be calculated some other way dependent on error calculations
    if (Lx2 / Nx2 < Lx1 / Nx1){
        dt = 1/4 * pow((Lx2 / Nx2),2);
    }
    if (dt == 0){
        dt = DBL_MIN;
        cout << "Uses smallest possible double, but it is still to big to garante for the error." << endl;
    }
    startX1 = Nx1 / 4;
    startX2 = Nx2 / 2;
    startX3 = Nx3 / 2;
    
    V = new double [Nx1 * Nx2 * Nx3];
    psi_r1 = new double [Nx1 * Nx2 * Nx3];
    psi_i1 = new double [Nx1 * Nx2 * Nx3];
    psi_r2 = new double [Nx1 * Nx2 * Nx3];
    psi_i2 = new double [Nx1 * Nx2 * Nx3];
    
    setV();
    makeInitState();
    finiteDifference();
}

void Schrodinger::contSim(string filename, unsigned int numOfTimesteps){
}

Schrodinger::Schrodinger(){
    V = nullptr;
    psi_r1 = nullptr;
    psi_i1 = nullptr;
    psi_r2 = nullptr;
    psi_i2 = nullptr;
}

Schrodinger::~Schrodinger(){
    delete [] V;
    delete [] psi_r1;
    delete [] psi_i1;
    delete [] psi_r2;
    delete [] psi_i2;
}

//PRIVATE MEMBER FUNCTIONS
void Schrodinger::setV(){
    V = new double [Nx1 * Nx2 * Nx3];
    switch (potential) {
        case FREE:
            setVtoZero();
            break;
        case CONST_BARRIER_1D:
            setVtoZero();
            for (int x1 = (Nx1/2 - VThickness*dx1); x1 < Nx1/2 + VThickness; x1++){
                V[x1] = V0;
            }
            break;
        case CONST_BARRIER_2D:
            setVtoZero();
            for (int x1 = (Nx1/2 - VThickness*dx1); x1 < Nx1/2 + VThickness; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                        V[Nx1*x2+x1] = V0;
                }
            }
            
            break;
        default:
            break;
    }
}

void Schrodinger::setVtoZero(){
    for (int x1 = 0; x1 < Nx1; x1++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x3 = 0; x3 < Nx3; x3++){
                V[Nx1*Nx2*x3+Nx1*x2+x1] = 0;
            }
        }
    }
}

void Schrodinger::makeInitState(){
    switch (probDistrb) {
        case GAUSSIAN_1D:
            for (int x1 = 0; x1 < Nx1; x1++){
                psi_r1[x1] = exp(-pow(dx1 * x1 - startX1, 2) / (2 * SDx1));
                psi_r2[x1] = 0;
                psi_i1[x1] = exp(-pow(dx1 * x1 - startX1, 2) / (2 * SDx1));
                psi_i2[x1] = 0;
            }
        case GAUSSIAN_2D:
            for (int x1 = 0; x1 < Nx1; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                    psi_r1[x1 * Nx2 + x2] = exp(-pow(dx1 * x1 - startX1, 2) / (2 * SDx1) - pow(dx1 * x2 - startX2, 2) / (2 * SDx2)) * cos(p * (dx1 * x1));
                    psi_r2[x1 * Nx2 + x2] = 0;
                    psi_i1[x1 * Nx2 + x2] = exp(-pow(dx1 * x1 - startX1, 2) / (2 * SDx1) - pow(dx2 * x2 - startX2, 2) / (2 * SDx2)) * sin(p * (dx1 * x1));
                    psi_i2[x1 * Nx2 + x2] = 0;
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
    FILE* plotFile = fopen((filename + "_plot").c_str(), "wb");
    double c1 = dt * hbar * hbar / 2 / m / dx1 / dx1;
    double c2 = dt / hbar;
    for (int t = 0; t < Nt; t++){
        for (int x = 1; x < Nx1 - 1; x++){
            psi_r2[x] = psi_r1[x] + (2 * c1 + c2 * V[x]) * psi_i1[x] - c1 * psi_i1[x + 1] - c1 * psi_i1[x - 1];
            psi_i2[x] = psi_i1[x] - (2 * c1 + c2 * V[x]) * psi_i1[x] + c1 * psi_i1[x + 1] + c1 * psi_i1[x - 1];
            psi_r1[x] = psi_r2[x] + (2 * c1 + c2 * V[x]) * psi_i2[x] - c1 * psi_i2[x + 1] - c1 * psi_i2[x - 1];
            psi_i1[x] = psi_i2[x] - (2 * c1 + c2 * V[x]) * psi_i2[x] + c1 * psi_i2[x + 1] + c1 * psi_i2[x - 1];
        }
        if (t % plotDensityT == 0){
            for (int x = 0; x < Nx1; x += plotDensityX1){
                if (numOfDim == 1){
                    fwrite(&psi_r1[x], sizeof(double), 1, plotFile);
                    fwrite(&psi_i1[x], sizeof(double), 1, plotFile);
                }
                double possibility = psi_r1[x] * psi_r1[x] + psi_i1[x] * psi_i1[x];
                fwrite(&possibility, sizeof(double), 1, plotFile);
            }
        }
    }
    fclose(plotFile);
    FILE* finalStateFile = fopen((filename + "_finalState").c_str(), "wb");
    fwrite(&psi_r1[0], sizeof(double), Nx1, finalStateFile);
    fwrite(&psi_i1[0], sizeof(double), Nx1, finalStateFile);
    fclose(finalStateFile);
    FILE* potentialFile = fopen((filename + "_potential").c_str(), "wb");
    fwrite(&V[0], sizeof(double), Nx1, potentialFile);
    fclose(potentialFile);
    writeVariablesToFile();
}

void Schrodinger::finiteDifference2D(){
    
}

void Schrodinger::finiteDifference3D(){
    
}

void Schrodinger::writeVariablesToFile(){
    ofstream finalStateFile(filename + "_variables.txt");
    finalStateFile << numOfDim << endl << Lx1 << endl << Lx2 << endl << Lx3 << endl << Nx1 << endl << Nx2 << endl << Nx3 << endl << Nt << endl << dx1 << endl << dx2 << endl << dx3 << endl << dt << endl << m << endl << p << endl << startX1 << endl << startX2 << endl << startX3 << endl << V0 << endl << VThickness << endl << situation << endl << potential << endl << probDistrb << endl <<  SDx1 << endl << SDx2 << endl << SDx3 << endl << plotDensityX1 << endl << plotDensityX2 << endl << plotDensityX3 << endl << plotDensityT << endl;
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
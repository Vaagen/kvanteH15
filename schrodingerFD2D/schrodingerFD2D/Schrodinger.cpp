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
    Nx1 = 1200;
    Nx2 = 1000;
    Nx3 = 1000;
    Nt = 100000;
    V0 = 0.0;
    VThickness = 0.0;
    m = 1;
    p = 1;
    
    SDx1 = 2 * Lx1 * Lx1 * Lx1;
    SDx2 = Lx2 * Lx2 * Lx2;
    SDx3 = Lx3 * Lx3 * Lx3;
    
    plotDensityX1 = 1;
    plotDensityX2 = 1;
    plotDensityX3 = 1;
    plotDensityT = Nt/500;
    
    // set the situation
    // this is where you add new situations (along with adding i new situation in the enum Situation)
    switch (situation) {
        case FREE_ELECTRON_1D:
            numOfDim = 1;
            m = pow(10, -30);
            potential = FREE;
            probDistrb = GAUSSIAN_1D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 2000 * m;
            Nt = 100000;
            break;
        case FREE_ELECTRON_2D:
            numOfDim = 2;
            m = pow(10, -30);
            potential = FREE;
            probDistrb = GAUSSIAN_2D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 10;
            break;
        case ELECTRON_CONST_BARRIER_1D:
            numOfDim = 1;
            m = pow(10, -30);
            potential = CONST_BARRIER_1D;
            V0 = pow(10, -30);
            VThickness = 0.0001;
            probDistrb = GAUSSIAN_1D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 10;
            Nt = 100000;
            break;
        case ELECTRON_TRIANGLE_1D:
            numOfDim = 1;
            m = pow(10, -30);
            potential = TRIANGLE_1D;
            V0 = pow(10, -30);
            VThickness = 0.0001;
            probDistrb = GAUSSIAN_1D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 10;
            Nt = 100000;
            // VdistanceToMax is distance from Lx1/2 to max of V, value is set in setV() under case: TRIANGLE_1D
            break;
        case ELECTRON_CONST_BARRIER_2D:
            numOfDim = 2;
            m = pow(10, -30);
            potential = CONST_BARRIER_2D;
            probDistrb = GAUSSIAN_2D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 10;
            break;
        default:
            break;
    }
    
    if (numOfDim == 1){
        Nx2 = 1;
        Lx2 = 1;
        Nx3 = 1;
        Lx3 = 1;
    } else if (numOfDim == 2){
        Nx3 = 1;
        Lx3 = 1;
    }
    dx1 = Lx1 / Nx1;
    dx2 = Lx2 / Nx2;
    dx3 = Lx3 / Nx3;
    k = 2 * 3.1415926535897 * 30 / Lx1; //p/hbar;
    startX1 = Nx1 / 4;
    startX2 = Nx2 / 2;
    startX3 = Nx3 / 2;
    
    V = new double [Nx1 * Nx2 * Nx3];
    psi_r1 = new double [Nx1 * Nx2 * Nx3];
    psi_i1 = new double [Nx1 * Nx2 * Nx3];
    psi_r2 = new double [Nx1 * Nx2 * Nx3];
    psi_i2 = new double [Nx1 * Nx2 * Nx3];
    
    setV();
    
    dt = hbar/(hbar * hbar/(2*m*dx1*dx1) + Vmax) * 0.001;//2.0 * pow(Lx1 / Nx1, 2); // should be calculated some other way dependent on error calculations
    if (dt == 0){
        dt = DBL_MIN;
        cout << "Uses smallest possible double as timestep, but it is still to big to garante for the error." << endl;
    }
    
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
            Vmax = 0;
            break;
        case CONST_BARRIER_1D:
            setVtoZero();
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                V[x1] = V0;
            }
            Vmax = V0;
            break;
        {case TRIANGLE_1D:
            setVtoZero();
            double VdistanceToMax = VThickness*0.5; // See Situation()
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VdistanceToMax/dx1; x1++){
                V[x1] = V0/(VdistanceToMax/dx1)*(x1-Nx1/2);
            }
            for (int x1 = (Nx1/2) + VdistanceToMax/dx1; x1 < Nx1/2 + VThickness/dx1; x1++){
                V[x1] = V0 - V0/((VThickness-VdistanceToMax)/dx1)* (x1-Nx1/2 - VdistanceToMax/dx1);
            }
            Vmax = V0;
            break;}
        case CONST_BARRIER_2D:
            setVtoZero();
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                        V[Nx1*x2+x1] = V0;
                }
            }
            Vmax = V0;
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
                psi_r1[x1] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1)) * cos(k * dx1 * x1);
                psi_r2[x1] = 0;
                psi_i1[x1] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1)) * sin(k * dx1 * x1);
                psi_i2[x1] = 0;
            }
            break;
        case GAUSSIAN_2D:
            for (int x1 = 0; x1 < Nx1; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                    psi_r1[x1 * Nx2 + x2] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1) - pow(dx1 * (x2 - startX2), 2) / (2 * SDx2)) * cos(p * (dx1 * x1));
                    psi_r2[x1 * Nx2 + x2] = 0;
                    psi_i1[x1 * Nx2 + x2] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1) - pow(dx2 * (x2 - startX2), 2) / (2 * SDx2)) * sin(p * (dx1 * x1));
                    psi_i2[x1 * Nx2 + x2] = 0;
                }
            }
            break;
            
        default:
            break;
    }
    normalizePsi();
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
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), "wb");
    FILE* plotPsiRFile = fopen((filename + "_plot_psi_r").c_str(), "wb");
    FILE* plotPsiIFile = fopen((filename + "_plot_psi_i").c_str(), "wb");
    double c1 = (dt / 2 / m / dx1 / dx1) * hbar;
    double c2 = dt / hbar;
    cout << c1 << endl;
    if (c1 == 0){
        cout << "c1 is 0 and the wavefunction will not change. You will probably need to make dt bigger." << endl;
    }
    for (int t = 0; t < Nt; t++){
        if (t % plotDensityT == 0){;
            for (int x = 0; x < Nx1; x += plotDensityX1){
                if (numOfDim == 1){
                    fwrite(&psi_r1[x], sizeof(double), 1, plotPsiRFile);
                    fwrite(&psi_i1[x], sizeof(double), 1, plotPsiIFile);
                }
                double possibility = psi_r1[x] * psi_r1[x] + psi_i1[x] * psi_i1[x];
                fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
            }
        }
        for (int x = 1; x < Nx1 - 1; x++){
            psi_r2[x] = psi_r1[x] + (2 * c1 + c2 * V[x]) * psi_i1[x] - c1 * psi_i1[x + 1] - c1 * psi_i1[x - 1];
            psi_i2[x] = psi_i1[x] - (2 * c1 + c2 * V[x]) * psi_r1[x] + c1 * psi_r1[x + 1] + c1 * psi_r1[x - 1];
            psi_r1[x] = psi_r2[x] + (2 * c1 + c2 * V[x]) * psi_i2[x] - c1 * psi_i2[x + 1] - c1 * psi_i2[x - 1];
            psi_i1[x] = psi_i2[x] - (2 * c1 + c2 * V[x]) * psi_r2[x] + c1 * psi_r2[x + 1] + c1 * psi_r2[x - 1];
        }
    }
    fclose(plotProbabilityFile);
    fclose(plotPsiRFile);
    fclose(plotPsiIFile);
    finalStore();
}

void Schrodinger::finiteDifference2D(){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), "wb");
    
    
    fclose(plotProbabilityFile);
    finalStore();
}

void Schrodinger::finiteDifference3D(){
    finalStore();
}

void Schrodinger::finalStore(){
    FILE* finalStateFile = fopen((filename + "_finalState").c_str(), "wb");
    fwrite(&psi_r1[0], sizeof(double), Nx1, finalStateFile);
    fwrite(&psi_i1[0], sizeof(double), Nx1, finalStateFile);
    fclose(finalStateFile);
    FILE* potentialFile = fopen((filename + "_potential").c_str(), "wb");
    fwrite(&V[0], sizeof(double), Nx1, potentialFile);
    fclose(potentialFile);
    probability = findProbability();
    writeVariablesToFile();
}

void Schrodinger::writeVariablesToFile(){
    ofstream finalStateFile(filename + "_variables.txt");
    finalStateFile << numOfDim << endl << Lx1 << endl << Lx2 << endl << Lx3 << endl << Nx1 << endl << Nx2 << endl << Nx3 << endl << Nt << endl << dx1 << endl << dx2 << endl << dx3 << endl << dt << endl << m << endl << p << endl << k << endl << startX1 << endl << startX2 << endl << startX3 << endl << V0 << endl << VThickness << endl << Vmax << endl << startEnergy << endl << finalEnergy << endl << finalProb << endl << situation << endl << potential << endl << probDistrb << endl <<  SDx1 << endl << SDx2 << endl << SDx3 << endl << plotDensityX1 << endl << plotDensityX2 << endl << plotDensityX3 << endl << plotDensityT << endl;
}

void Schrodinger::normalizePsi(){
    for (int x1 = 0; x1 < Nx1; x1++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x3 = 0; x3 < Nx3; x3++){
                psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] /= findProbability();
                psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] /= findProbability();
            }
        }
    }
}

double Schrodinger::findProbability(){
    double probability;
    for (int x1 = 0; x1 < Nx1; x1++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x3 = 0; x3 < Nx3; x3++){
                probability += (psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] + psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_i1[Nx1*Nx2*x3+Nx1*x2+x1]) * dx1 * dx2 * dx3;
            }
        }
    }
    return probability;
}

double findEnergy(){
    return 1.0;
}
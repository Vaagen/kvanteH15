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
    Lx2 = 0.0005;
    Lx3 = 0.001;
    Nx1 = 1200;
    Nx2 = 500;
    Nx3 = 1000;
    Nt = 1000;
    V0 = 0.0;
    VThickness = 0.0;
    m = 1;
    p = 1;
    startEnergy = 0.0;
    finalEnergy = 0.0;
    
    SDx1 = 2 * Lx1 * Lx1 * Lx1;
    SDx2 = Lx2 * Lx2 * Lx2;
    SDx3 = Lx3 * Lx3 * Lx3;
    
    plotSpacingX1 = 1;
    plotSpacingX2 = 1;
    plotSpacingX3 = 1;
    
    double numOfFrames = 500;
    
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
            Nx1 = 1000;
            m = pow(10, -30);
            potential = FREE;
            probDistrb = GAUSSIAN_2D;
            plotSpacingX1 = 10;
            plotSpacingX2 = 5;
            numOfFrames = 10;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            Nt = 10000;
            p = 10;
            break;
        case ELECTRON_CONST_BARRIER_1D:
            numOfDim = 1;
            m = pow(10, -30);
            potential = CONST_BARRIER_1D;
            V0 = pow(10, -28);
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
        case ELECTRON_MULTIPLE_SLIT_2D:
            numOfDim = 2;
            m = pow(10, -30);
            potential = MULTIPLE_SLIT_2D;
            probDistrb = GAUSSIAN_2D;
            //SDx1 = SDx1;
            //SDx2 = SDx2;
            p = 10;
            // The following variables are set in setV() under case: MULTIPLE_SLIT_2D
            // slitNumber, number of slits in barrier
            // slitWidth, width of each slit
            // slitDistance, distance between each slit
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
    plotSpacingT = Nt/numOfFrames;
    
    V = new double [Nx1 * Nx2 * Nx3];
    psi_r1 = new double [Nx1 * Nx2 * Nx3];
    psi_i1 = new double [Nx1 * Nx2 * Nx3];
    psi_r2 = new double [Nx1 * Nx2 * Nx3];
    psi_i2 = new double [Nx1 * Nx2 * Nx3];
    
    setV();
    
    dt = hbar/(hbar * hbar/(2*m*dx1*dx1) + Vmax) * 0.001;//2.0 * pow(Lx1 / Nx1, 2); // should be calculated some other way dependent on error calculations
    if (numOfDim ==2){
        dt = 0.0001;
    }
    if (dt == 0){
        dt = DBL_MIN;
        cout << "Uses smallest possible double as timestep, but it is still to big to garante for the error." << endl;
    }
    
    makeInitState();
    startEnergy = findEnergy();
    
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
        {case MULTIPLE_SLIT_2D:
            int slitNumber = 2;
            double slitWidth = Lx2/15;
            double slitDistance = Lx2/20;
            setVtoZero();
            // Making constant potential barrier
            for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                for (int x2 = 0; x2 < Nx2; x2++){
                    V[Nx1*x2+x1] = V0;
                }
            }
            // Making slits in constant potential barrier
            int slitsPlaced = 0;
            double nextSlitX2 = 0; // to keep track of where to place next slit
            // Placing first slit, if odd number of slits.
            if (slitNumber % 2 == 0){
                nextSlitX2 = slitDistance/(dx2*2);
            }else if (slitsPlaced % 2 == 1){
                for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                    for (int x2 = Nx2/2 - slitWidth/(dx2*2); x2 < Nx2/2 + slitWidth/(dx2*2); x2++){
                        V[Nx1*x2+x1] = 0;
                    }
                }
                slitsPlaced++;
                nextSlitX2 = slitWidth/2 + slitDistance;
            }
            // Placing remaining slits, two at a time
            while (slitsPlaced < slitNumber && (nextSlitX2 + slitWidth < Lx2/2)){
                for (int x1 = (Nx1/2); x1 < Nx1/2 + VThickness/dx1; x1++){
                    for (int x2 = nextSlitX2/dx2; x2 < nextSlitX2/dx2 + slitWidth/dx2; x2++){
                        V[Nx1*(Nx2/2 + x2) +x1] = 0;
                        V[Nx1*(Nx2/2 - x2) +x1] = 0;
                    }
                }
                slitsPlaced += 2;
                nextSlitX2 += slitWidth + slitDistance;
            }
            Vmax = V0;
            break;}
        default:
            break;
    }
}

void Schrodinger::setVtoZero(){
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
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
            for (int x2 = 0; x2 < Nx2; x2++){
                for (int x1 = 0; x1 < Nx1; x1++){
                    psi_r1[Nx1*x2 + x1] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1) - pow(dx2 * (x2 - startX2), 2) / (2 * SDx2)) * cos(p * (dx1 * x1));
                    psi_r2[Nx1*x2 + x1] = 0;
                    psi_i1[Nx1*x2 + x1] = exp(-pow(dx1 * (x1 - startX1), 2) / (2 * SDx1) - pow(dx2 * (x2 - startX2), 2) / (2 * SDx2)) * sin(p * (dx1 * x1));
                    psi_i2[Nx1*x2 + x1] = 0;
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
    finalStore();
    cout << "The final energy is " << finalEnergy/startEnergy << " times the start energy." << endl;
    cout << "The final probability of finding the particle is: " << finalProb << endl;
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
        double possibility;
        if (t % plotSpacingT == 0){;
            for (int x = 0; x < Nx1; x += plotSpacingX1){
                fwrite(&psi_r1[x], sizeof(double), 1, plotPsiRFile);
                fwrite(&psi_i1[x], sizeof(double), 1, plotPsiIFile);
                possibility = psi_r1[x] * psi_r1[x] + psi_i1[x] * psi_i1[x];
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
}

void Schrodinger::finiteDifference2D(){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), "wb");
    double c1x1 = hbar * dt / dx1 / dx1;
    double c1x2 = hbar * dt / dx2 / dx2;
    c1x1 = 0.008;
    c1x2 = 0.008;
    double c2 = dt / hbar;
    cout << "c1x1: " << c1x1 << endl;
    cout << "c1x2: " << c1x2 << endl;
    if (c1x1 == 0 || c1x2 == 0){
        cout << "c1x1 or c1x2 is 0 and the wavefunction will not change. You will probably need to make dt bigger." << endl;
    }
    for (int t = 0; t < Nt; t++){
        double possibility;
        if (t % plotSpacingT == 0){
            for (int x2 = 1; x2 < Nx2 - 1; x2+=plotSpacingX2){
                for (int x1 = 1; x1 < Nx1 - 1; x1+=plotSpacingX1){
                    possibility = psi_r1[x2*Nx1 + x1] * psi_r1[x2*Nx1 + x1] + psi_i1[x2*Nx1 + x1] * psi_i1[x2*Nx1 + x1];
                    fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
                }
            }
        }
        int i;
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*x2+x1;
                psi_r2[i] = psi_r1[i] + (2*c1x1 + 2*c1x2 + c2*V[i])*psi_i1[i] - c1x1*(psi_i1[i+1] + psi_i1[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]);
                psi_i2[i] = psi_i1[i] - (2*c1x1 + 2*c1x2 + c2*V[i])*psi_r1[i] + c1x1*(psi_r1[i+1] + psi_r1[i-1]) + c1x2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]);
                psi_r1[i] = psi_r2[i] + (2*c1x1 + 2*c1x2 + c2*V[i])*psi_i2[i] - c1x1*(psi_i2[i+1] + psi_i2[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i2[i-Nx1]);
                psi_i1[i] = psi_i2[i] - (2*c1x1 + 2*c1x2 + c2*V[i])*psi_r2[i] + c1x1*(psi_r2[i+1] + psi_r2[i-1]) + c1x2*(psi_r2[i+Nx1] + psi_r2[i-Nx1]);
            }
        }
    }
    fclose(plotProbabilityFile);
    finalStore();
}

void Schrodinger::finiteDifference3D(){
    FILE* plotProbabilityFile = fopen((filename + "_plot_probability").c_str(), "wb");
    double c1x1 = hbar * dt / dx1 / dx1;
    double c1x2 = hbar * dt / dx2 / dx2;
    double c1x3 = hbar * dt / dx3 / dx3;
    double c2 = dt / hbar;
    for (int t = 0; t < Nt; t++){
        double possibility;
        if (t % plotSpacingT == 0){
            for (int x3 = 0; x3 < Nx3; x3+= plotSpacingX3){
                for (int x2 = 0; x2 < Nx2; x2+=plotSpacingX2){
                    for (int x1 = 0; x1 < Nx1; x1+=plotSpacingX1){
                        possibility = psi_r1[Nx1*Nx2*x3 + Nx1*x2 + x1] * psi_r1[Nx1*Nx2*x3 + Nx1*x2 + x1] + psi_i1[Nx1*Nx2*x3 + Nx1*x2 + x1] * psi_i1[Nx1*Nx2*x3 + Nx1*x2 + x1];
                        fwrite(&possibility, sizeof(double), 1, plotProbabilityFile);
                    }
                }
            }
        }
        int i;
        for (int x3 = 1; x3 < Nx3 - 1; x3++){
            for (int x2 = 1; x2 < Nx2 - 1; x2++){
                for (int x1 = 1; x1 < Nx1 - 1; x1++){
                    i = Nx1*Nx2*x3+Nx1*x2+x1;
                    psi_r2[i] = psi_r1[i] + (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_i1[i] - c1x1*(psi_i1[i+1] + psi_i1[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i1[i-Nx1]) - c1x3*(psi_i1[i+Nx1*Nx2] + psi_i1[i-Nx1*Nx2]);
                    psi_i2[i] = psi_i1[i] - (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_r1[i] + c1x1*(psi_r1[i+1] + psi_r1[i-1]) + c1x2*(psi_r1[i+Nx1] + psi_r1[i-Nx1]) + c1x3*(psi_r1[i+Nx1*Nx2] + psi_r1[i-Nx1*Nx2]);
                    psi_r1[i] = psi_r2[i] + (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_i2[i] - c1x1*(psi_i2[i+1] + psi_i2[i-1]) - c1x2*(psi_i1[i+Nx1] + psi_i2[i-Nx1]) - c1x3*(psi_i2[i+Nx1*Nx2] + psi_i2[i-Nx1*Nx2]);
                    psi_i1[i] = psi_i2[i] - (2*c1x1 + 2*c1x2 + 2*c1x3 + c2*V[i])*psi_r2[i] + c1x1*(psi_r2[i+1] + psi_r2[i-1]) + c1x2*(psi_r2[i+Nx1] + psi_r2[i-Nx1]) + c1x3*(psi_r2[i+Nx1*Nx2] + psi_r2[i-Nx1*Nx2]);
                }
            }
        }
    }
    fclose(plotProbabilityFile);
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
    finalProb = findProbability();
    finalEnergy = findEnergy();
    writeVariablesToFile();
}

void Schrodinger::writeVariablesToFile(){
    ofstream finalStateFile(filename + "_variables.txt");
    finalStateFile << numOfDim << endl << Lx1 << endl << Lx2 << endl << Lx3 << endl << Nx1 << endl << Nx2 << endl << Nx3 << endl << Nt << endl << dx1 << endl << dx2 << endl << dx3 << endl << dt << endl << m << endl << p << endl << k << endl << startX1 << endl << startX2 << endl << startX3 << endl << V0 << endl << VThickness << endl << Vmax << endl << startEnergy << endl << finalEnergy << endl << finalProb <<  endl << situation << endl << potential << endl << probDistrb << endl <<  SDx1 << endl << SDx2 << endl << SDx3 << endl << plotSpacingX1 << endl << plotSpacingX2 << endl << plotSpacingX3 << endl << plotSpacingT << endl;

}

void Schrodinger::normalizePsi(){
    double probability = findProbability();
    //cout << probability << endl;
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] /= sqrt(probability);
                psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] /= sqrt(probability);
            }
        }
    }
    
}

double Schrodinger::findProbability(){
    double probability = 0.0;
    for (int x3 = 0; x3 < Nx3; x3++){
        for (int x2 = 0; x2 < Nx2; x2++){
            for (int x1 = 0; x1 < Nx1; x1++){
                probability += (psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_r1[Nx1*Nx2*x3+Nx1*x2+x1] + psi_i1[Nx1*Nx2*x3+Nx1*x2+x1] * psi_i1[Nx1*Nx2*x3+Nx1*x2+x1]) * dx1 * dx2 * dx3;
            }
        }
    }
    return probability;
}

double Schrodinger::findEnergy(){
    if (numOfDim == 1) {
        return findEnergy1D();
    } else if (numOfDim == 2){
        return findEnergy2D();
    } else {
        return findEnergy3D();
    }
}

double Schrodinger::findEnergy1D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double energy = 0.0;
    int i;
    for (int x1 = 1; x1 < Nx1 - 1; x1++){
        i = x1;
        energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + (V[i]-2*cx1)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1;
    }
    return energy;
}

double Schrodinger::findEnergy2D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double cx2 = -hbar * hbar / 2 / m / dx2 / dx2;
    double energy = 0.0;
    int i;
    for (int x2 = 1; x2 < Nx2 - 1; x2++){
        for (int x1 = 1; x1 < Nx1 - 1; x1++){
            i = Nx1*x2 + x1;
            energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + cx2*(psi_r1[i]*(psi_r1[i+Nx1]+psi_r1[i-Nx1]) + psi_i1[i]*(psi_i1[i+Nx1]+psi_i1[i-Nx1])) + (V[i]-2*cx1-2*cx2)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1;
        }
    }
    return energy;
}

double Schrodinger::findEnergy3D(){
    double cx1 = -hbar * hbar / 2 / m / dx1 / dx1;
    double cx2 = -hbar * hbar / 2 / m / dx2 / dx2;
    double cx3 = -hbar * hbar / 2 / m / dx3 / dx3;
    double energy = 0.0;
    int i;
    for (int x3 = 1; x3 < Nx3 - 1; x3++){
        for (int x2 = 1; x2 < Nx2 - 1; x2++){
            for (int x1 = 1; x1 < Nx1 - 1; x1++){
                i = Nx1*Nx2*x3 + Nx1*x2 + x1;
                energy += (cx1*(psi_r1[i]*(psi_r1[i+1]+psi_r1[i-1]) + psi_i1[i]*(psi_i1[i+1]+psi_i1[i-1])) + cx2*(psi_r1[i]*(psi_r1[i+Nx1]+psi_r1[i-Nx1]) + psi_i1[i]*(psi_i1[i+Nx1]+psi_i1[i-Nx1]))  + cx3*(psi_r1[i]*(psi_r1[i+Nx1*Nx2]+psi_r1[i-Nx1*Nx2]) + psi_i1[i]*(psi_i1[i+Nx1*Nx2]+psi_i1[i-Nx1*Nx2])) + (V[i]-2*cx1-2*cx2-2*cx3)*(psi_r1[i]*psi_r1[i] + psi_i1[i]*psi_i1[i]))*dx1;
            }
        }
    }
    return energy;
}
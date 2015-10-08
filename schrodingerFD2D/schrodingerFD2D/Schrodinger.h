//
//  Schrodinger.h
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#ifndef __schrodingerFD__Schrodinger__
#define __schrodingerFD__Schrodinger__

#include <stdio.h>
#include <cmath>
#include <string>
#include <iostream>
#include <vector>

using namespace std;

enum Situation{FREE_ELECTRON_1D, FREE_ELECTRON_2D, ELECTRON_CONST_BARRIER_1D, ELECTRON_CONST_BARRIER_2D};
enum Potential{FREE,CONST_BARRIER_1D,CONST_BARRIER_2D};
enum ProbDistrb{GAUSSIAN_1D, GAUSSIAN_2D}; //probability distribution for initial position


// uses finite difference method to simulate the schrodinger equation in 2 dimensions from a variety initial conditil
// If you wish to add new situations or tweek the parameters, make a new constant in the enum "Situation" (and situationString) and a new switch-case in the private member function "runSituation" and make your changes inside this new case. New potentials can be make the same way with "Potentials" and "setV".
class Schrodinger{
public:
    // run the simulation with 'situation' and store it under filename
    void run(Situation situation, string filename);
    
    // continue previous simulation
    void contSim(string filename, unsigned int numOfTimesteps);
    
    // sets the referanse potential
    void setV0(double V0){this->V0 = V0;}
    
    // sets the referanse width
    void setVThickness(double VThickness){this->VThickness = VThickness;}
    
    Schrodinger();
    ~Schrodinger();
private:
// VARIABLES FOR SIMULATION
    const double hbar = 1; //1.0545718 * pow(10, -34);
    
    string filename;
    unsigned int numOfDim;
    double Lx1;
    double Lx2;
    double Lx3;
    unsigned int Nx1;
    unsigned int Nx2;
    unsigned int Nx3;
    unsigned int Nt; //number of timesteps
    double dx1;
    double dx2;
    double dx3;
    double dt;
    double m;
    double p;
    double k;
    double startX1;
    double startX2;
    double startX3;
    double V0; // referance potential
    double VThickness; // referance thickness
    
    Situation situation;
    Potential potential;
    ProbDistrb probDistrb;
    
    double SDx1; // uncertity in position
    double SDx2; // uncertity in position
    double SDx3; // uncertity in position
    
    unsigned int plotDensityX1; //steps between x1 plotted
    unsigned int plotDensityX2; //steps between x2 plotted
    unsigned int plotDensityX3; //steps between x3 plotted
    unsigned int plotDensityT; //steps between t plotted
    
    double* V;
    double* psi_r1; //real part of wave function
    double* psi_i1; //imagenary part of wave funciton
    double* psi_r2; //real part of wave function
    double* psi_i2; //imagenary part of wave funciton

    
// MEMBER FUNCTIONS
    // spesifies the potential
    void setV();
    // sets V to zero everywhere
    void setVtoZero();
    
    // setting initial state
    void makeInitState();
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference();
    
    void finiteDifference1D();
    
    void finiteDifference2D();
    
    void finiteDifference3D();
    
    void writeVariablesToFile();
};




#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */

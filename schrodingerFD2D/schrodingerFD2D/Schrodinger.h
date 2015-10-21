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

enum Situation{FREE_ELECTRON_1D, FREE_ELECTRON_2D, ELECTRON_CONST_BARRIER_1D, ELECTRON_TRIANGLE_1D, ELECTRON_CONST_BARRIER_2D, ELECTRON_MULTIPLE_SLIT_2D, ELECTRON_BALL_2D};
enum Potential{FREE,CONST_BARRIER_1D,TRIANGLE_1D, CONST_BARRIER_2D, MULTIPLE_SLIT_2D, CIRCLE_2D, BALL_2D};
enum ProbDistrb{GAUSSIAN_1D, GAUSSIAN_2D}; //probability distribution for initial position


// uses finite difference method to simulate the schrodinger equation in 2 dimensions from a variety initial conditions
// If you wish to add new situations or tweek the parameters, make a new constant in the enum "Situation" (and situationString) and a new switch-case in the private member function "runSituation" and make your changes inside this new case. New potentials can be make the same way with "Potentials" and "setV".
class Schrodinger{
public:
    // run the simulation with 'situation' and store it under filename
    void run(Situation situation, string filename);
    
    // continue previous simulation
    void continueSimulation(string filename, unsigned int numOfTimesteps);
    
    // sets the referanse potential
    void setV0(double V0){this->V0 = V0;}
    
    // sets the referanse width
    void setVThickness(double VThickness){this->VThickness = VThickness;}
    
    Schrodinger();
    ~Schrodinger();
private:
// VARIABLES FOR SIMULATION
    const double hbar = 1.0e0; // 1.0545718 * pow(10, -34);
    
    string filename;
    unsigned int numOfDim;
    double Lx1;
    double Lx2;
    double Lx3;
    int Nx1;
    int Nx2;
    int Nx3;
    int Nt; //number of iterations (timesteps = 3*Nt)
    double dx1;
    double dx2;
    double dx3;
    double dt;
    double m;
    double p;
    double k;
    int startX1;
    int startX2;
    int startX3;
    double V0; // referance potential
    double VThickness; // potential thickness
    double Vmax;
    double startEnergy;
    double finalEnergy;
    double finalProb;
    
    Situation situation;
    Potential potential;
    ProbDistrb probDistrb;
    
    double SDx1; // uncertity in position
    double SDx2; // uncertity in position
    double SDx3; // uncertity in position
    
    int plotSpacingX1; //steps between x1 plotted
    int plotSpacingX2; //steps between x2 plotted
    int plotSpacingX3; //steps between x3 plotted
    int plotSpacingT; //steps between t plotted
    
    int numOfFrames; // Nt / plotSpacingT
    
    double* V;
    double* psi_r1; //real part of wave function
    double* psi_i1; //imagenary part of wave funciton
    double* psi_r2; //real part of wave function
    double* psi_i2; //imagenary part of wave funciton
    double* psi_r3; //real part of wave function
    double* psi_i3; //imagenary part of wave funciton
    
// MEMBER FUNCTIONS
    // spesifies the potential
    void setV();
        // sets V to zero everywhere
        void setVtoZero();
    
    // setting initial state
    void makeInitState();
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference(bool newSimulation);
        void finiteDifference1D(char* fileOpenType);
        void finiteDifference2D(char* fileOpenType);
        void finiteDifference3D(char* fileOpenType);
            void finalStore();
    
    void writeVariablesToFile();
    
    void normalizePsi();
    double findProbability();
    double findEnergy();
        double findEnergy1D();
        double findEnergy2D();
        double findEnergy3D();
    
    void loadVaiables();
    void loadFinalState();
};




#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */

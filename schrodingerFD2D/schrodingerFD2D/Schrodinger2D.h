//
//  Schrodinger2D.h
//  schrodingerFD2D
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#ifndef __schrodingerFD2D__Schrodinger2D__
#define __schrodingerFD2D__Schrodinger2D__

#include <stdio.h>
#include <cmath>
#include <string>

using namespace std;

enum Potential{FREE};
enum InitState{};


// uses finite difference method to simulate the schrodinger equation in 2 dimensions from a variety initial conditil
class Schrodinger2D{
public:
    
    // run the simulation
    // let you choose from a list of initial contitions
    void run();
    
    // make a inital condition
    // let you choose type of particle, speed, probability distribusjon and potential
    void makeInitCondition();
    
    // continue previous simulation
    void contSim(string filename, unsigned int numOfTimesteps);
    
private:
    // VARIABLES FOR SIMULATION
    double Lx;
    double Ly;
    unsigned int Nx;
    unsigned int Ny;
    double dt;
    unsigned int Nt; //number of timesteps
    double t;
    const double hbar = 1.0545718 * pow(10, -34);
    double m;
    double* V;
    double* psi_r; //real part of wave function
    double* psi_i; //imagenary part of wave funciton
    unsigned int plotDensity; //steps between plots
    
    // MEMBER FUNCTIONS
    // spesifies the potential
    void setV(Potential type);
    
    // setting initial state
    void initState(string type);
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference();
};




#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */

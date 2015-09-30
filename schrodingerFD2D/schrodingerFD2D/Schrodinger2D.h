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
#include <iostream>
#include <vector>

using namespace std;

enum Potential{FREE};
enum Situation{FREE_ELECTRON};
vector<string> situationString{"free_electron"};
enum ProbDistrb{GAUSSIAN}; //probability distribution for initial position


// uses finite difference method to simulate the schrodinger equation in 2 dimensions from a variety initial conditil
// If you wish to add new situations or tweek the parameters, make a new constant in the enum "Situation" (and situationString) and a new switch-case in the private member function "runSituation" and make your changes inside this new case. New potentials can be make the same way with "Potentials" and "setV".
class Schrodinger2D{
public:
    // run the simulation
    // let you choose from a list of situations
    void run();
    
    // run the simulation with 'situation'
    void runSituation(Situation situation);
    
    // continue previous simulation
    void contSim(string filename, unsigned int numOfTimesteps);
    
private:
// VARIABLES FOR SIMULATION
    const double hbar = 1.0545718 * pow(10, -34);

    double Lx;
    double Ly;
    unsigned int Nx;
    unsigned int Ny;
    unsigned int Nt; //number of timesteps
    double dx;
    double dy;
    double dt;
    double m;
    double p;
    double startX;
    double startY;
    
    Potential potential;
    Situation situation;
    ProbDistrb probDistrb;
    double SDx; // uncertity in position
    double SDy; // uncertity in position
    
    unsigned int plotDensityX; //steps between x plotted
    unsigned int plotDensityY; //steps between y plotted
    unsigned int plotDensityT; //steps between t plotted
    
    double* V;
    double* psi_r; //real part of wave function
    double* psi_i; //imagenary part of wave funciton

    
// MEMBER FUNCTIONS
    // spesifies the potential
    void setV();
    
    // setting initial state
    void makeInitState();
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference();
};




#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */

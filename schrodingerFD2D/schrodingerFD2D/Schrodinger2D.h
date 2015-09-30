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

class Schrodinger2D{
public:
    // run the simulation
    // let you choose type of particle, speed, probability distribusjon and potential
    void run();
    // continue simulation
    void continue
private:
    // VARIABLES FOR SIMULATION
    double Lx;
    double Ly;
    int Nx;
    int Ny;
    double dt;
    int Nt;
    const double hbar = 1.0545718 * pow(10, -34);
    double** psi_r; //real part of wave function
    double** psi_i; //imagenary part of wave funciton
    
    // MEMBER FUNCTIONS
    //
    // spesifies the potential
    void V(string type);
    
    // setting initial state
    void initState(string type);
    
    // calculating the time evolution and storing it in a tekst file
    void finiteDifference();
};

#endif /* defined(__schrodingerFD2D__Schrodinger2D__) */

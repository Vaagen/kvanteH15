//
//  main.cpp
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include <iostream>
#include "Schrodinger.h"
#include <ctime>
#include <cstdio>
using namespace std;

int main(int argc, const char * argv[]) {
    clock_t t = clock() / CLOCKS_PER_SEC;
    Schrodinger s;
    //s.run(FREE_ELECTRON_1D, "test_free_electron");
    //s.run(ELECTRON_CONST_BARRIER_1D, "test_free_electron");
    s.run(FREE_ELECTRON_2D, "test_free_electron");
    //s.run(ELECTRON_TRIANGLE_1D, "test_free_electron");
    //s.run(ELECTRON_MULTIPLE_SLIT_2D, "test_free_electron");
    //s.continueSimulation("test_free_electron", 1000);
    long time = static_cast<long>(clock() - t) / CLOCKS_PER_SEC;
    cout << "The simulation used " << time / 60 << " minuttes and " << time % 60 << " seconds." << endl;
}
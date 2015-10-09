//
//  main.cpp
//  schrodingerFD
//
//  Created by Håkon Austlid Taskén on 30.09.15.
//  Copyright (c) 2015 Håkon Austlid Taskén. All rights reserved.
//

#include <iostream>
#include "Schrodinger.h"
using namespace std;

int main(int argc, const char * argv[]) {    Schrodinger s;
    //s.run(FREE_ELECTRON_1D, "test_free_electron");
    s.run(ELECTRON_CONST_BARRIER_1D, "test_free_electron");
    //s.run(ELECTRON_TRIANGLE_1D, "test_free_electron_");
    
}
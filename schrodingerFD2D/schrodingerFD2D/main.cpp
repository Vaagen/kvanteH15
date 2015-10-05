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

int main(int argc, const char * argv[]) {
    int a = 1000;
    int b = 1000;
    int c = 1000;
    double* test = new double [a*b*c];
    for (int i = 0; i < a; i++){
        for (int j = 0; j < b; j++){
            for (int k = 0; k < c; k++){
                test[c*b*i + c*j + k] = i*j;
            }
        }
    }
    delete [] test;
    cout << "yes" << endl;
}
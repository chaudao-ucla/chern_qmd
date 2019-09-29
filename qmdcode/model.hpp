#ifndef model_hpp
#define model_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>

#include "vec3.hpp"
#include "util.hpp"

class Neighbor {
public:
    int idx;
    vec3 delta;
    double dr;
};

const int Max_N_neighbor = 28; // 6 + 12 + 8

class Super_Cell {
public:
    int idx;
    int x, y, z;
    double r_min[3];
    double r_max[3];
    vector<int> atoms;
    
    int N_neighbor;
    
    Super_Cell *nn[Max_N_neighbor];
    
//    Super_Cell *nn1[6];
//    Super_Cell *nn2[12];
//    Super_Cell *nn3[8];
};

struct AtomPair {
    int index1;
    int index2;
    vec3 delta; // r_2 - r_1
};

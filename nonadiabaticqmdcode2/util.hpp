//
//  util.hpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#ifndef util_hpp
#define util_hpp

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <random>

//#define ARMA_DONT_USE_WRAPPER
#include <armadillo>

//#include "vec3.hpp"

#include <vector>

// Thin wrapper that adds bounds checking to std::vector
template<typename T>
class Vec : public std::vector<T> {
public:
    // using std::vector<T>::vector;
    // Be explicit to avoid annoying Xcode errors
    Vec<T>(): std::vector<T>() {}
    Vec<T>(int size): std::vector<T>(size) {}
    Vec<T>(int size, T const& val): std::vector<T>(size, val) {}
    Vec<T>(std::initializer_list<T> l): std::vector<T>(l) {}
    template <class InputIterator>
    Vec<T>(InputIterator first, InputIterator last): std::vector<T>(first, last) {}
    
    T& operator[](int i) {
        return std::vector<T>::at(i);
    }
    
    const T& operator[](int i) const {
        return std::vector<T>::at(i);
    }
};



const double PI = acos(-1.0);

typedef float flt;
// typedef double flt;
typedef std::complex<flt> cx_flt;
typedef std::complex<double> cx_double;

typedef std::mt19937 RNG;


const cx_double _I(0, 1);

void select_sort(double *a, int _L);

void select_sort(double *a, int _L, int *idx);

void init_srand(void);

inline double rand1(void) {
    //return ((double) rand())/((double) RAND_MAX);
    return drand48();
}

int mod(int x, int m);

double fermi_density(double x, double kT, double mu);





#endif /* util_hpp */

//
//  test.cpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include "test.hpp"
//#include "model.hpp"
#include "util.hpp"

using namespace std;

void test_all(void) {
    
    test_sparse_matrix();
    
    exit(1);
}

void test_sparse_matrix(void) {
    int dim = 40000;
    arma::sp_cx_mat A(dim, dim);
    arma::sp_cx_mat B(dim, dim);
    arma::sp_cx_mat C(dim, dim);
    
    std::random_device seed;
    
    RNG rng = RNG(seed());
    
    std::uniform_real_distribution<double> rd;
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    
    for(int i=0; i<50000; i++) {
        int m = (int) (dim * rd(rng));
        int n = (int) (dim * rd(rng));
        A(m, n) = rn(rng);
    }
    
    
    for(int i=0; i<dim-1; i++) B(i, i+1) = 0.23 * i * (i-1);
    
    C = A * B;
    
    cout << C << endl;
}


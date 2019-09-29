//
//  test.cpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <armadillo>

#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"
#include "test.hpp"

using namespace std;

void test_all(void) {
    
    test_sparse_matrix();
    
    exit(1);
}

void test_Bott_idx(void) {
    
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.00000001;                // Langevin damping coefficient
    int N_atoms = 100;
    double kT = 0.001;
    double density = 4.0 * 0.03;
    
    double l_box = pow(N_atoms / density, 1./3.);
    
    cout << "l_box = " << l_box << endl;
    
    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT, gamma);
    atoms.set_BD_type(1);
    
    cout << "TI-TB parameters: " << endl;
    double M = 2.5;
    double lambda = 0.5;
    double t2 = 0.25;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    atoms.init_potential();
    
    double Bott_idx = atoms.compute_Bott_index();
    
    cout << "Bott index = " << Bott_idx << endl;

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


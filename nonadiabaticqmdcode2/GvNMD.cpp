//
//  main.cpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <iostream>
#include "util.hpp"
#include "model.hpp"
#include "test.hpp"

using namespace std;

int main(int argc, const char * argv[]) {
    // insert code here...
    std::cout << "Hello, World!\n";
    
    double U_i              = argc > 1 ? atof(argv[1]) : 0.01;
    double U_f              = argc > 2 ? atof(argv[2]) : 5.0;
    double rs               = argc > 3 ? atof(argv[3]) : 0.85;

    
    double kT_start         = argc > 4 ? atof(argv[4]) : 0.1;
    double kT_stop          = argc > 5 ? atof(argv[5]) : 0.005;
    int n_annealing         = argc > 6 ? atoi(argv[6]) : 25;
    double r_annealing      = argc > 7 ? atof(argv[7]) : 0.96;
    
    double dlt_t            = argc > 8 ? atof(argv[8]) : 0.01;
    int max_steps           = argc > 9 ? atoi(argv[9]) : 200000;
    
    
    int N_atoms  = 10;

    double l_box = pow((4. * M_PI / 3.) * ((double) N_atoms), 1./3.) * rs;     // simulation box size

    cout << "U_i = " << U_i << "\t U_f = " << U_f << endl;
    cout << "dt = " << dlt_t << "\t max_steps = " << max_steps << endl;
    
    cout << "rs = " << rs << endl;
    cout << "l_box = " << l_box << endl;

    double r_max = 7.5;

    init_srand();
    //test_all();
    
    TB_model system(rs, N_atoms, r_max);
    
    cout << "L = " << system.L << endl;
    
    system.GA_kT_start = kT_start;
    system.GA_kT_stop = kT_stop;
    system.GA_n_annealing = n_annealing;
    system.GA_r_annealing = r_annealing;
    system.filling = 0.5;
    
    system.h0 = 24.;
    system.alpha = 2.0;
    system.v0 = 100;
    system.beta = 3.53;
    system.mass = 1;
    
    system.r_annealing = 0.97;

    /*
    system.kT = 0.001;
    system.gamma = 0.00000001;
    system.QMD_simulation(10000, 0.01, 0.000000000001);
    exit(1);
    */
    
    //test_all();

    /*
    system.init_random();
    system.build_cell_atoms_list();
    system.compute_neighbor_lists();
    
    system.init_bare_hopping();
    
    system.U = U_i;
    
    system.compute_GA_solution();
    */
    
    
    //system. plot_binding_curve();
    
    system.kT_init_annealing = 0.5;
    system.kT_end_annealing = 0.005;

    system.simulate_quench(max_steps, dlt_t, U_i, U_f);
    
    return 0;
}

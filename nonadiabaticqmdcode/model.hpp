//
//  model.hpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

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

using namespace std;

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

class TB_model {
public:
    
    double rs, L;
    int Ns;
    int dim;
    
    vec3 system_box;
    
    double h0, alpha;
    double v0, beta;
    double r_max;
    
    double U;
    
    double filling;
    double mu;
    double kT;
    
    double mass;
    
    arma::mat Mpp, Npp, Upp;
    
    std::mt19937 rng;
    std::random_device seed;
    
    static constexpr int N_nn1 = 4;
    class Atom {
    public:
        int idx;
        
        vec3 position;
        
        vec3 velocity;
        vec3 force;
        
        int cell_idx;
        
        vector<Neighbor> neighbor_list;
        
        cx_double R;
        cx_double Delta;
        
        double lambda;
        
        Atom(void) {
            neighbor_list.clear();
        };
        
    } *atom;
    
    arma::cx_mat Hopping;
    
    arma::sp_cx_mat Hamiltonian;
    arma::cx_mat Density_Mat;
    arma::cx_mat Phi;
    
    arma::mat position;
    arma::mat momentum;
    
    int N_c;
    int L_c;
    double c_size;
    Super_Cell *cell;
    
    
    TB_model(double r0, int N0, double r_cutoff) {
        
        rs = r0;
        Ns = N0;
        
        r_max = r_cutoff;
        
        L = pow(4. * M_PI / 3., 1./3.) * rs * pow((double) Ns, 1./3);
        system_box.x = L;
        system_box.y = L;
        system_box.z = L;
        
        cout << "linear system size L = " << L << endl;
        
        dim = Ns;
        
        atom = new Atom[Ns];
        
        Hopping = arma::cx_mat(dim, dim);
        
        Hamiltonian = arma::sp_cx_mat(dim, dim);
        Density_Mat = arma::cx_mat(dim, dim);
        
        Phi = arma::cx_mat(3, dim);
        
        position = arma::mat(3, dim);
        momentum = arma::mat(3, dim);
        
        init_pp_matrices();
        
        build_super_cells(0.5 * r_max);
        
//        unsigned int seed2 = (unsigned int) (rand1() * 100000000000.);
//        rng = RNG(seed() + seed2);
        
        rng = RNG(seed());
        
        langevin_dym = 1;
    };
    
    ~TB_model(void) {
        
        delete [] atom;
        delete [] cell;
    };
    
    // ==================================================
    
    int index_cell(int x, int y, int z) {
        return L_c * L_c * z + L_c * y + x;
    };
    
    void build_super_cells(double r_c);
    void build_cell_atoms_list(void);
    void compute_neighbor_lists(void);

    vector<AtomPair> find_all_pairs(double cutoff);
    void update_atoms(arma::mat &);
    
    // ==================================================

    double hopping_func(double r) {
        return -h0 * exp(-alpha * r);
    }
    
    vec3 dt_dr(vec3 delta) {
        double r = delta.norm();
        
        return h0 * alpha * exp(-alpha * r) * delta / r;
    }
    
    
    double phi(double r) {
        
        return v0 * exp(-beta * r);
    }
    
    double dphi_dr(double r) {
        return -v0 * beta * exp(-beta * r);
    }
    
    vec3 dphi_dr(vec3 delta) {
        double r = delta.norm();
        
        return -v0 * beta * exp(-beta * r) * delta / r;
    }

    // ==================================================

    void init_random(void);
    void init_from_file(void);
    
    double wrapDelta1d(double dx, double boxLength);
    vec3 wrapDelta(vec3 delta, vec3 boxLengths);
    double wrapPosition1d(double x, double bdsLo, double bdsHi);
    vec3 wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi);
    void move_atoms_pbc(void);

    
    void compute_hopping_matrix(void);
    void compute_hopping_matrix_full(void);

    
    // ==================================================

    void init_pp_matrices(void);
    
    
    void reset_GA_parameters(void) {
        for(int i=0; i<Ns; i++) {
            atom[i].R = 1.0;
            atom[i].lambda = U;
        }
    };
    
    void clear_hopping(void) {
        for(int i=0; i<Ns; i++) atom[i].neighbor_list.clear();
    };
    
    // ==================================================
    
    arma::sp_cx_mat build_Hamiltonian(void);
    arma::sp_cx_mat build_Hamiltonian(arma::mat &);
    
    void compute_fermi_level(arma::vec & eigE);
    arma::cx_mat compute_density_matrix(arma::sp_cx_mat &);
    arma::cx_mat compute_density_matrix(arma::cx_mat &);
    arma::cx_mat compute_density_matrix(arma::mat &);
    
    arma::cx_mat compute_tV_density_matrix(arma::sp_cx_mat &);
    void compute_GA_solution(void);
    
    void compute_Deltas(arma::cx_mat &);
    void compute_Deltas(arma::cx_mat &, arma::mat &);
    
    void compute_renormalizations(arma::cx_mat &, arma::cx_mat &);
    arma::cx_mat solve_Phi(arma::cx_mat &);
    double adjust_lambda(arma::cx_mat &, arma::cx_mat &);
    
    arma::cx_mat init_uncorrelated_Phi(arma::cx_mat &);
    
    void statistics_GA(int ndata, double U);
    
    
    // ========================================
    
    double GA_kT_start;
    double GA_kT_stop;
    int GA_n_annealing;
    double GA_r_annealing;
    
    double avg_R = 0, sgm_R = 0;
    double avg_nd = 0, sgm_nd = 0;
    double avg_np = 0, sgm_np = 0;
    double avg_d = 0, sgm_d = 0;
    
    double E_e = 0;
    double E_p = 0;
    double E_k = 0;
    
    void compute_electronic_energy(void);
    void compute_pair_potential(void);
    void compute_kinetic_energy(void);
    
    vec3 V_cm;
    void compute_cm_velocity(void);
    void adjust_cm_velocity(void);
    
    void analyze_data(void);
    
    void test(void);

    // ==========================================
        
    arma::mat compute_dX(arma::mat &);
    arma::mat compute_dP(arma::mat &, arma::cx_mat &);
    
    arma::cx_mat compute_dF(arma::cx_mat &, arma::cx_mat &);
    void integrate_EOM_RK4(double dt);
    void integrate_EOM_RK4_v2(double dt);
    
    void init_config(void);
    void simulate_quench(int max_steps, double dt, double U_i, double U_f);
    void simulate_noneq(int max_steps, double dt, double U_i, double U_f);

    
    // ==========================================

    
    double gamma;
    int langevin_dym = 1;
    
    void move_atoms(double dt);
    void compute_forces(arma::cx_mat & Dm);
    void compute_forces_full(arma::cx_mat & Dm);
    void integrate_forces(double dt);
    void integrate_Langevin(double dt);
    void step_NVT(double dt);
    
    double compute_E_elec();
    double compute_E_pair();
    double compute_E_kin();
    
    double compute_E_pair_full();
    
    double r_annealing;
    
    void plot_binding_curve(void);
    
    double kT_init_annealing;
    double kT_end_annealing;
    
    void QMD_simulation(int N_steps, double dt);
    
};

#endif /* model_hpp */

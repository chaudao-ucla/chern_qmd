//
//  tbmd.hpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef tbmd_hpp
#define tbmd_hpp

#include <armadillo>
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

#include "potential.hpp"

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

class Simple_TB_system 
{
    public:
        int numAtoms;

        int N_c;
        int L_c;
        double c_size;
        Super_Cell *cell;

        vector<vec3> position;
        vector<vec3> velocity;
        vector<vec3> force;

        vector<vec3> pos_init;
        vector<vec3> v_init;

        vec3 cm_velocity;

        double mass;
        double kT;
        double gamma; // Langevin damping

        double filling_fraction; // used in equilibrium case
        int num_filling;
        double mu;

        double kT_electron;

        arma::vec eigE; // current eigen-energies of the Hamiltonian
        arma::cx_mat eig_mode; // current eigenstates of the hamiltonian
        arma::cx_mat proj_states; // U_m(t) for non-equilibrium case
        arma::cx_vec fermi_energies; // Original fermi-energies

        // pointer to tight-binding potential models:
        // ToyModelSOribtal *potential;
        TI_2D_atom *potential;
        
        std::mt19937 rng;
        std::random_device seed;
        
        int dimH;
        
        //fkpm::EnergyScale energyScale;
        arma::cx_mat Hamiltonian;
        arma::cx_mat rho;
        
        arma::cx_mat U_mat, W_mat;
        
        arma::cx_mat Vx_mat, Vy_mat;
        
        vec3 system_box;
        double volume;
        
        // constructor:
        Simple_TB_system(int nAtoms, double atomMass, double filling, double temperature, double gm);
        
        // initialization:
        void init_potential(void);
        
        // ============================================================================
        // initialization of various lattices:
        void init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);
        void init_from_file(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);
        void init_from_file_vel(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi);
        
        // ============================================================================
        // set parameters:
        int boundary_type;   // 0 : Open BC,  1 : Periodic BC
        int langevin_dym;
        
        void set_BD_type(int t) {boundary_type = t;};   // default = 0
        void set_LD_dynm(int t) {langevin_dym = t;};    // default = 1
        
        // ============================================================================
        vector<AtomPair> pairs;

        void find_all_pairs(void);
        void find_all_pairs(double cutoff);
        vector<AtomPair> get_pairs(vector<vec3> pos);

        // ==================================================
    
        int index_cell(int x, int y, int z) 
        {
            return L_c * L_c * z + L_c * y + x;
        };

        void build_super_cells(double r_c);
        void build_cell_atoms_list(void);
        void compute_neighbor_lists(void);
    
        // ==================================================
        
        // ============================================================================
        void built_Hamiltonian(void);
        
        //fkpm::EnergyScale computeEnergyScale(double extra, double tolerance);
        
        // ============================================================================
        // for single-particle density matrix:
        double compute_chemical_potential(void);
        double compute_avg_density(double x);
        arma::cx_mat compute_density_matrix(arma::cx_mat & H_elec);
        void compute_density_matrix(void);
        void proj_to_rho(arma::cx_mat const & proj, arma::cx_mat & denmat);

        arma::cx_mat evolve_density_matrix_RK1(arma::cx_mat & H_elec, 
                arma::cx_mat & proj_states, arma::cx_vec & fermi_energies, double dt);

        void evolve_density_matrix_RK1(double dt);
        
        // ============================================================================
        // MD-related routines
        void move_atoms(double dt);
        void compute_forces(void);
        void integrate_forces(double dt);
        void integrate_Langevin(double dt);
        void step_NVT(double dt);

        void evolve_NVT_RK1(double dt);
        void evolve_NVT_RK4(double dt);
        void evolve_NVT_RK4_v2(double dt);
        void move_atoms_RK1(double dt);
        void integrate_forces_RK1(double dt);

        void scale(double factor);
        //void compute_forces_elec(void);
        
        // ============================================================================
        // compute energies
        double ee, ep, ek;
        
        double e_gap;
        
        double e_kin(void);
        double e_elec(void);
        double e_pair(void);
        
        double compute_Eg(void);
        
        // measurement
        double virial;
        
        
        // compute Bott index:
        double compute_Bott_index(void);
        
        // compute Hall conductivity
        void compute_velocity_matrix(void);
        double compute_sigma_xy(void);
        
        void save_configuration(string const filename);
        
        void plot_eigen_mode(string const filename, int m);
        void print_eig_E(string const filename);
        
        void test_1(void);
};

#endif /* tbmd_hpp */

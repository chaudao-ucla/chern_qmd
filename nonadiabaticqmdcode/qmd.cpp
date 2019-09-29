//
//  qmd.cpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/4/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include <vector>
#include "model.hpp"
#include "analysis.hpp"

using namespace std;

vector<AtomPair> TB_model::find_all_pairs(double cutoff) {
  //  (vector<vec3> const& pos, double cutoff, vec3 boxLengths) {

    vector<AtomPair> pairs;
    for(int i=0; i<Ns; i++) {
        for(int j=i+1; j<Ns; j++) {
            vec3 delta = wrapDelta(atom[j].position - atom[i].position, system_box);
            if(delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {i, j, delta});
        }
    }
    
    return pairs;
}

void TB_model::move_atoms(double dt) {
    for(int i=0; i<Ns; i++) {
        vec3 dlt = 0.5 * (atom[i].force / mass) * dt;
        dlt += atom[i].velocity;
        atom[i].velocity = dlt;  // velecity at t + 0.5*dt
        dlt *= dt;
        atom[i].position += dlt;
    }
}

void TB_model::integrate_forces(double dt) {
    for(int i=0; i<Ns; i++) {
        atom[i].velocity += 0.5 * (atom[i].force / mass) * dt;
    }
}

void TB_model::integrate_Langevin(double dt) {
    //RNG rng(seed());
    std::normal_distribution<double> rd;    // default mean = 0, var = 1
    
    double alpha2 = exp(-gamma * dt);
    double sigma2 = sqrt((1-pow(alpha2, 2))*kT/mass);
    for(int i=0; i<Ns; i++) {
        atom[i].velocity = alpha2 * atom[i].velocity + sigma2 * vec3(rd(rng), rd(rng), rd(rng));
    }
}

void TB_model::compute_forces(arma::cx_mat & Dm) {
    for(int i=0; i<Ns; i++) {
        
        atom[i].force = vec3(0, 0, 0);
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;
            vec3 delta = (*k).delta;  // wrapDelta(atom[j].position - atom[i].position, system_box);
            
            atom[i].force += dphi_dr(delta) + 4. * real(Dm(j, i)) * dt_dr(delta);
            
        }
    }

}

void TB_model::compute_forces_full(arma::cx_mat & Dm) {
    for(int i=0; i<Ns; i++) {
        
        atom[i].force = vec3(0, 0, 0);
        
        for(int j=0; j<Ns; j++) {

            if(j != i) {
            
                vec3 delta = wrapDelta(atom[j].position - atom[i].position, system_box);
                atom[i].force += dphi_dr(delta) + 4. * real(Dm(j, i)) * dt_dr(delta);
            
            }
        }
    }
    
}

void TB_model::step_NVT(double dt) {
    
    // velocity verlet integration with Langevin damping
    
    move_atoms(dt);
    move_atoms_pbc();
    
    build_cell_atoms_list();
    compute_neighbor_lists();
    compute_hopping_matrix();
    
//    compute_hopping_matrix_full();
    
    Density_Mat = compute_density_matrix(Hopping);
    
    compute_forces(Density_Mat);
//    compute_forces_full(Density_Mat);

    integrate_forces(dt);
    if(langevin_dym == 1) integrate_Langevin(dt);
}

double TB_model::compute_E_elec(void) {
    
    cx_double sum = 0;
    
    /*
    for(int i=0; i<Ns; i++) {
        cx_double tmp;
        
        tmp = Hopping(i, i);
        sum += (tmp) * Density_Mat(i, i);
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            int j = (*k).idx;
            tmp = Hopping(i, j);
            sum += tmp * Density_Mat(j, i);
        }
        
    }
    */
    
    auto rr = Hopping * Density_Mat;
    sum = trace(rr).real();

   
    E_e = 2. * sum.real() / ((double) Ns);

    return E_e;
}

double TB_model::compute_E_pair(void) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            //            int j = (*k).idx;
            double dr = (*k).dr;
            
            sum += phi(dr);
        }
    }
    
    E_p = 0.5 * sum / ((double) Ns);
    
    return E_p;
}

double TB_model::compute_E_pair_full(void) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
       
        for(int j=0; j<Ns; j++) {
            
            if(j != i) {
        
                vec3 delta = wrapDelta(atom[j].position - atom[i].position, system_box);
                double dr = delta.norm();
                sum += phi(dr);
            
            }
        }
    }
    
    E_p = 0.5 * sum / ((double) Ns);
    
    return E_p;
}

double TB_model::compute_E_kin(void) {
    double sum = 0;
    for(int i=0; i<Ns; i++) sum += atom[i].velocity.norm2();
    E_k = (0.5 * mass * sum) / ((double) Ns);
    
    return E_k;
}

void TB_model::QMD_simulation(int N_steps, double dt) {
    
    MD_Data mddata;
    MD_Data accu_data;
    int n_therm_save = 50;
    int nx = 0;
    
    double r_min = 0.0002;
    double r_max = 0.5 * L;
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    double D_max = 0.3;
    int num_hs_D = 180;
    double dlt_D = D_max/((double) num_hs_D);
    mddata.init_hist_D(num_hs_D, dlt_D);
    
    double R_max = 1.0;
    int num_hs_R = 150;
    double dlt_R = R_max/((double) num_hs_R);
    mddata.init_hist_R(num_hs_R, dlt_R);

    
    ofstream fs("p.dat");
    
    for(int i=0; i<Ns; i++) {
        atom[i].velocity = vec3(0, 0, 0);
    }

    init_random();
    //init_from_file();
    
    build_cell_atoms_list();
    compute_neighbor_lists();
    compute_hopping_matrix();
    Density_Mat = compute_density_matrix(Hopping);
    compute_forces(Density_Mat);
    
//    compute_hopping_matrix_full();
//    Density_Mat = compute_density_matrix(Hopping);
//    compute_forces_full(Density_Mat);
    
    kT = kT_init_annealing;
    
    for(int r=0; r < N_steps; r++) {
        
        if(r % 100 == 0)
            cout << "r = " << r << "\t kT = " << kT << endl;
        
        kT = kT_end_annealing + (kT_init_annealing - kT_end_annealing) * pow(1. - (((double) r) / ((double) N_steps)), 4.);
        
        step_NVT(dt);
        
        compute_E_elec();
        compute_E_pair();
//        compute_E_pair_full();
        compute_E_kin();
        
        fs << r * dt << '\t' << E_e << '\t' << E_p << '\t' << E_k << '\t';
        fs << endl;
        
        if(r % n_therm_save == 0) {
            
            mddata.basic_measurements(*this);
            
            auto pairs = find_all_pairs(1.5 * r_max);
            
            mddata.compute_g_r(pairs, Ns, L*L*L, r_min, r_max);
            
            if(r > 1000) {
                average_MD_Data(accu_data, mddata, nx);
                nx++;
            }
            accu_data.print_g_r("fr.dat");
            
        }
    }
    
    fs.close();

}

void TB_model::plot_binding_curve(void) {
    
    ofstream fs("eb.dat");
    double dr = 0.001;
    for(double r = 0.001; r<5; r+= dr) {
        fs << hopping_func(r) << '\t' << phi(r) << '\t' << -2 * hopping_func(r) + phi(r) << endl;
    }
    fs.close();
}

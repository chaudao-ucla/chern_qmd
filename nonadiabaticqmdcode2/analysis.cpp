//
//  analysis.cpp
//  GA_solver-TBMD
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "model.hpp"
#include "analysis.hpp"

void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc) {
    
    if(nc == 0) {
        cout << "resetting accumulated MD_data." << endl;
        accu_data.reset();
        accu_data.init_g_r(mddata.num_r, mddata.dr);
        accu_data.init_hist_D(mddata.num_hD, mddata.dlt_D);
        accu_data.init_hist_R(mddata.num_hR, mddata.dlt_R);
    }
    
    accu_data.M1_ee = (accu_data.M1_ee * nc + mddata.M1_ee) / (nc + 1.);
    accu_data.M2_ee = (accu_data.M2_ee * nc + mddata.M2_ee) / (nc + 1.);
    
    accu_data.M1_ep = (accu_data.M1_ep * nc + mddata.M1_ep) / (nc + 1.);
    accu_data.M2_ep = (accu_data.M2_ep * nc + mddata.M2_ep) / (nc + 1.);
    
    accu_data.M1_ek = (accu_data.M1_ek * nc + mddata.M1_ek) / (nc + 1.);
    accu_data.M2_ek = (accu_data.M2_ek * nc + mddata.M2_ek) / (nc + 1.);

    
    for(unsigned m=0; m<accu_data.g_r.size(); m++) {
        accu_data.g_r[m] = (accu_data.g_r[m] * nc + mddata.g_r[m]) / (nc + 1.);
    }
    
    for(unsigned m=0; m<accu_data.hs_D.size(); m++) {
        accu_data.hs_D[m] = (accu_data.hs_D[m] * nc + mddata.hs_D[m]) / (nc + 1.);
    }
    
    for(unsigned m=0; m<accu_data.hs_R.size(); m++) {
        accu_data.hs_R[m] = (accu_data.hs_R[m] * nc + mddata.hs_R[m]) / (nc + 1.);
    }
    
}

void MD_Data::basic_measurements(TB_model & system) {
    
    system.compute_electronic_energy();
    system.compute_pair_potential();
    system.compute_kinetic_energy();
    
    double ee = system.E_e;
    double ep = system.E_p;
    double ek = system.E_k;
    
    M1_ee = ee;
    M2_ee = pow(ee, 2);
    
    M1_ep = ep;
    M2_ep = pow(ep, 2);
    
    M1_ek = ek;
    M2_ek = pow(ek, 2);
        
}

void MD_Data::compute_g_r(vector<AtomPair> pairs, int Nm, double volume, double r_min, double r_max) {
    
    for(int m=0; m<num_r; m++) g_r[m] = 0;
    
    for(int i=0; i<pairs.size(); i++) {
        double rij = pairs[i].delta.norm();
        if(rij < r_max && rij > r_min) {
            int m = (int) ((rij - r_min)/dr);
            
            if(m<num_r && m>=0)
                g_r[m]++;
        }
    }
    
    for(int m=0; m<num_r; m++) {
        double rm = (m + 0.5) * dr;
        g_r[m] *= volume / (2. * 3.1415926 * dr * pow(Nm * rm, 2));
    }
}

void MD_Data::compute_hist_double_occupancy(TB_model & system, double d_max) {
    double d_d = d_max / ((double) num_hD);
    
    for(int m=0; m<num_hD; m++) hs_D[m] = 0;
    int N_tot = 0;
    
    for(int i=0; i<system.Ns; i++) {
        double dc = real(conj(system.Phi(2, i)) * system.Phi(2, i));
        int m = (int) (dc/d_d);
        
        if(m>=0 && m<num_hD) {
            N_tot++;
            hs_D[m]++;
        }
    }
    
    for(int m=0; m<num_hD; m++) {
        hs_D[m] /= ((N_tot + 1.e-5) * d_d);
    }
}

void MD_Data::compute_hist_renormalization(TB_model & system, double r_max) {
    double d_r = r_max / ((double) num_hR);
    
    for(int m=0; m<num_hR; m++) hs_R[m] = 0;
    int N_tot = 0;
    
    for(int i=0; i<system.Ns; i++) {
        double Rm = real(system.atom[i].R * conj(system.atom[i].R));
        int m = (int) (Rm/d_r);
        
        if(m>=0 && m<num_hR) {
            N_tot++;
            hs_R[m]++;
        }
    }
    
    for(int m=0; m<num_hR; m++) {
        hs_R[m] /= (N_tot * d_r);
    }
}


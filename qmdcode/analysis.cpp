//
//  analysis.cpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"

void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc) {
    
    if(nc == 0) {
        cout << "resetting accumulated MD_data." << endl;
        accu_data.reset();
        accu_data.init_g_r(mddata.num_r, mddata.dr);
        
    }
    
    accu_data.M1_ee = (accu_data.M1_ee * nc + mddata.M1_ee) / (nc + 1.);
    accu_data.M2_ee = (accu_data.M2_ee * nc + mddata.M2_ee) / (nc + 1.);
    
    accu_data.M1_ep = (accu_data.M1_ep * nc + mddata.M1_ep) / (nc + 1.);
    accu_data.M2_ep = (accu_data.M2_ep * nc + mddata.M2_ep) / (nc + 1.);
    
    accu_data.M1_ek = (accu_data.M1_ek * nc + mddata.M1_ek) / (nc + 1.);
    accu_data.M2_ek = (accu_data.M2_ek * nc + mddata.M2_ek) / (nc + 1.);
    
    accu_data.E1 = (accu_data.E1 * nc + mddata.E1) / (nc + 1.);
    accu_data.E2 = (accu_data.E2 * nc + mddata.E2) / (nc + 1.);
    
    accu_data.M1_Eg = (accu_data.M1_Eg * nc + mddata.M1_Eg) / (nc + 1.);
    accu_data.M2_Eg = (accu_data.M2_Eg * nc + mddata.M2_Eg) / (nc + 1.);
    
    accu_data.M1_Bt = (accu_data.M1_Bt * nc + mddata.M1_Bt) / (nc + 1.);
    accu_data.M2_Bt = (accu_data.M2_Bt * nc + mddata.M2_Bt) / (nc + 1.);
    
    accu_data.M1_sgm = (accu_data.M1_sgm * nc + mddata.M1_sgm) / (nc + 1.);
    accu_data.M2_sgm = (accu_data.M2_sgm * nc + mddata.M2_sgm) / (nc + 1.);
    
    for(unsigned m=0; m<accu_data.g_r.size(); m++) {
        accu_data.g_r[m] = (accu_data.g_r[m] * nc + mddata.g_r[m]) / (nc + 1.);
    }
    
}

void MD_Data::basic_measurements(Simple_TB_system & system) {
    
    double ee = system.e_elec();
    double ep = system.e_pair();
    double ek = system.e_kin();
    
    
    M1_ee = ee;
    M2_ee = pow(ee, 2);
    
    M1_ep = ep;
    M2_ep = pow(ep, 2);
    
    M1_ek = ek;
    M2_ek = pow(ek, 2);
    
    E1 = ee + ep + ek;
    E2 = pow(E1, 2);
    
    double e_gap = system.compute_Eg();
    
    M1_Eg = e_gap;
    M2_Eg = pow(e_gap, 2);
    
    double bott_idx = system.compute_Bott_index();
    
    M1_Bt = bott_idx;
    M2_Bt = pow(bott_idx, 2);
    
    double sgm_xy = system.compute_sigma_xy();
    
    M1_sgm = sgm_xy;
    M2_sgm = pow(sgm_xy, 2);
    
}

void MD_Data::compute_g_r(vector<AtomPair> pairs, int Nm, double volume, double r_min, double r_max, int dim) {
    
    for(int m=0; m<num_r; m++) g_r[m] = 0;
    
    for(int i=0; i<pairs.size(); i++) {
        double rij = pairs[i].delta.norm();
        if(rij < r_max && rij > r_min) {
            int m = (int) ((rij - r_min)/dr);
            
            if(m<num_r && m>=0)
                g_r[m]++;
        }
    }
    
    if(dim == 3) {
        for(int m=0; m<num_r; m++) {
            double rm = (m + 0.5) * dr;
            g_r[m] *= volume / (2. * 3.1415926 * dr * pow(Nm * rm, 2));
        }
    } else if(dim == 2) {
        for(int m=0; m<num_r; m++) {
            double rm = (m + 0.5) * dr;
            g_r[m] *= volume / (3.1415926 * dr * rm * pow(Nm, 2));
        }
    }
}


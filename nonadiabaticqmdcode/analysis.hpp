//
//  analysis.hpp
//  GA_solver-TBMD
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#ifndef analysis_hpp
#define analysis_hpp

#include <stdio.h>
#include "util.hpp"
#include "model.hpp"

using namespace std;

class Measurement {
public:
    double dt;
    int numSteps;
    int numAtoms;
    int tr;
    
    int active;
    int enabled;
    
    Vec<vec3> ref_velocity;
    Vec<double> autocorr_v;
    Vec<double> accum_diff;
    double diffusion;
    
    Measurement(void) { active = 0; enabled = 0; };
    
    Measurement(double time_step, int N_dt, int N_atoms) {
        dt = time_step;
        numSteps = N_dt;    // must be an odd number !!!
        numAtoms = N_atoms;
        
        ref_velocity = Vec<vec3>(numAtoms);
        autocorr_v = Vec<double>(numSteps);
        accum_diff = Vec<double>(numSteps);
        
        active = 0;
        enabled = 0;
    };
    
    void start(TB_model & system) {
        
        active = 1;
        enabled = 1;
        tr = 0;
        
        double mass = system.mass;
        for(int i=0; i<system.Ns; i++) {
            ref_velocity[i].x = system.momentum(0, i) / mass;
            ref_velocity[i].y = system.momentum(1, i) / mass;
            ref_velocity[i].z = system.momentum(2, i) / mass;
        }
        
    };
    
    void record_data(TB_model & system) {
        
        double mass = system.mass;
        if(active == 1 && tr < numSteps) {
            
            // velocity autocorrelation measurement:
            double sum = 0;
            for(int i=0; i<system.Ns; i++) {
                sum += (system.momentum(0, i) * ref_velocity[i].x + system.momentum(1, i) * ref_velocity[i].y + system.momentum(2, i) * ref_velocity[i].z) / mass;
                
            }
            autocorr_v[tr] = sum / (3. * system.Ns);
            
            tr++;
        }
        
        if(tr == numSteps) {
            
            active = 0;
            
        }
        
    };
    
    void compute_diffusion_coef(void) {
        double sum = 0;
        for(int i=0; i<autocorr_v.size(); i++) {
            sum += autocorr_v[i];
            accum_diff[i] = sum * dt;
        }
        diffusion = sum * dt;
        
    };
    
};

class MD_Data {
public:
    
    double M1_ee, M2_ee;
    double M1_ep, M2_ep;
    double M1_ek, M2_ek;
    
    vector<double> g_r;
    int num_r;
    double dr;
    
    vector<double> hs_D;
    int num_hD;
    double dlt_D;
    
    vector<double> hs_R;
    int num_hR;
    double dlt_R;

    
    MD_Data(void) {
        reset();
    };
    
    ~MD_Data(void) {
        
    };
    
    void reset(void) {
        M1_ee = M2_ee = M1_ep = M2_ep = M1_ek = M2_ek = 0;
    };
    
    
    void init_g_r(int size, double dlt_r) {
        
        g_r.resize(size);
        num_r = size;
        dr = dlt_r;
        
        for(int i=0; i<size; i++) g_r[i] = 0;
    };
    
    void print_g_r(string const filename) {
        
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for(unsigned m=0; m<g_r.size(); m++) fs << (m + 0.5) * dr << '\t' <<  g_r[m] << endl;
        fs.close();
        
    };
    
    void init_hist_D(int size, double delta_D) {
        hs_D.resize(size);
        num_hD = size;
        dlt_D = delta_D;
        
        for(int i=0; i<size; i++) hs_D[i] = 0;
    };
    
    void init_hist_R(int size, double delta_R) {
        hs_R.resize(size);
        num_hR = size;
        dlt_R = delta_R;
        
        for(int i=0; i<size; i++) hs_R[i] = 0;
    };
    
    void print_hist_D(string const filename) {
        
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for(unsigned m=0; m<hs_D.size(); m++) fs << (m + 0.5) * dlt_D << '\t' <<  hs_D[m] << endl;
        fs.close();
        
    };
    
    void print_hist_R(string const filename) {
        
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(10);
        for(unsigned m=0; m<hs_R.size(); m++) fs << (m + 0.5) * dlt_R << '\t' <<  hs_R[m] << endl;
        fs.close();
        
    };
    
    void basic_measurements(TB_model & system);
    
    void compute_g_r(vector<AtomPair> pairs, int Nm, double volume, double r_min, double r_max);
    
    void compute_hist_double_occupancy(TB_model & system, double d_max);
    
    void compute_hist_renormalization(TB_model & system, double r_max);

    // ==============================================================
    // for velocity autocorrelation function and diffusion coefficient

    int nv;
    
    int numSteps;
    double dt;
    
    double M1_D, M2_D;
    Vec<double> autocorr_v;
    Vec<double> accum_diff;

    void init_time_series(int num_steps, double time_step) {
        
        nv = 0;
        
        numSteps = num_steps;
        dt = time_step;
        
        M1_D = M2_D = 0;
        
        autocorr_v = Vec<double>(num_steps);
        accum_diff = Vec<double>(num_steps);
        
        for(int k=0; k<num_steps; k++) {
            autocorr_v[k] = accum_diff[k] = 0;
        }
    };
    
    void average_ms_data(Measurement & ms) {
        
        M1_D = (M1_D * nv + ms.diffusion) / (nv + 1.);
        M2_D = (M2_D * nv + pow(ms.diffusion, 2)) / (nv + 1.);
        
        for(int m=0; m<ms.numSteps; m++) {
            
            autocorr_v[m] = (autocorr_v[m] * nv + ms.autocorr_v[m]) / (nv + 1.);
            accum_diff[m] = (accum_diff[m] * nv + ms.accum_diff[m]) / (nv + 1.);
            
        }
        
        nv++;
    };
    
    void print_auto_correlation(std::string const str) {
        
        std::ofstream fs;
        fs.open(str, std::ios::out);
        
        for(int m=0; m<numSteps; m++) {
            fs << m * dt << '\t';
            fs << autocorr_v[m] << '\t' << accum_diff[m] << '\t';
            fs << endl;
        }
        fs.close();
        
    };
    
    
    void print_data(std::string const str, double param) {
        
        std::ofstream fs;
        
        fs.open(str, std::ios::out);
        
        fs << param << '\t';
        fs << M1_ee << '\t' << sqrt(M2_ee - pow(M1_ee, 2)) << '\t';
        fs << M1_ep << '\t' << sqrt(M2_ep - pow(M1_ep, 2)) << '\t';
        fs << M1_ek << '\t' << sqrt(M2_ek - pow(M1_ek, 2)) << '\t';
        fs << M1_D << '\t'  << sqrt(M2_D - pow(M1_D, 2)) << '\t';
        fs << endl;
        
        fs.close();
        
    };

    
};


void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc);


#endif /* analysis_hpp */

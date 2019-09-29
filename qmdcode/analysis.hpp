//
//  analysis.hpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef analysis_hpp
#define analysis_hpp

#include <stdio.h>
#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"

using namespace std;

class MD_Data {
public:
    
    double M1_ee, M2_ee;
    double M1_ep, M2_ep;
    double M1_ek, M2_ek;
    
    double E1, E2;
    
    double M1_Eg, M2_Eg;
    
    double M1_Bt, M2_Bt;
    
    double M1_sgm, M2_sgm;
    
    vector<double> g_r;
    int num_r;
    double dr;
    
    
    int numSteps;
    double dt;
    
    
    MD_Data(void) {
        reset();
    };
    
    ~MD_Data(void) {
        
    };
    
    void reset(void) {
        M1_ee = M2_ee = M1_ep = M2_ep = M1_ek = M2_ek = 0;
        E1 = E2 = 0;
        M1_Eg = M2_Eg = 0;
        M1_Bt = M2_Bt = 0;
        M1_sgm = M2_sgm = 0;
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

    void print_g_r(std::ofstream & fs)
    {
        for(unsigned m=0; m<g_r.size(); m++) fs << (m + 0.5) * dr << '\t' <<  g_r[m] << endl;
    }
    
    void basic_measurements(Simple_TB_system & system);
    
    void compute_g_r(vector<AtomPair> pairs, int Nm, double volume, double r_min, double r_max, int dim);
    
    void print_data(string const filename, double param) {
        
        std::ofstream fs;
        
        fs.open(filename.c_str(), ios::out);
        fs.precision(12);
        
        fs << param << '\t';
        fs << M1_ee << '\t' << sqrt(fabs(M2_ee - pow(M1_ee, 2))) << '\t';
        fs << M1_ep << '\t' << sqrt(fabs(M2_ep - pow(M1_ep, 2))) << '\t';
        fs << M1_ek << '\t' << sqrt(fabs(M2_ek - pow(M1_ek, 2))) << '\t';
        fs << E1 << '\t' << (E2 - E1 * E1) << '\t';
        fs << M1_Eg << '\t' << sqrt(fabs(M2_Eg - pow(M1_Eg, 2))) << '\t';
        fs << M1_Bt << '\t' << sqrt(fabs(M2_Bt - pow(M1_Bt, 2))) << '\t';
        fs << M1_sgm << '\t' << sqrt(fabs(M2_sgm - pow(M1_sgm, 2))) << '\t';
        fs << endl;
        
        fs.close();
        
    };
    
};


void average_MD_Data(MD_Data & accu_data, MD_Data & mddata, int nc);

#endif /* analysis_hpp */

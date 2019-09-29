//
//  potential.hpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef potential_hpp
#define potential_hpp

#include <stdio.h>
#include <vector>
#include "util.hpp"
#include <armadillo>


using namespace std;

inline int positive_mod(int i, int n) {
    return (i%n + n) % n;
}

// Polynomial that splices together functions f and g between points x0 and x1, of the form
//     a0 + a1 (x - x0) + a2 (x - x0)^2 + a3 (x - x0)^3
// Coefficients a0..a3 are selected so that the spliced function and its derivatives are continuous.
struct SplicePolynomial {
    double a0, a1, a2, a3, x0;
    SplicePolynomial(double f0, double df0, double x0, double x1, double g1, double dg1);
    double eval(double x);
    double deriv(double x);
};

struct AtomPair {
    int index1;
    int index2;
    vec3 delta; // r_2 - r_1
};
std::ostream& operator<< (std::ostream& os, AtomPair const& p);



// Replace with fast neighbor finding algorithm
vec3 wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi);
vec3 wrapDelta(vec3 p, vec3 temp);
vector<AtomPair> allPairs(vector<vec3> const& pos, double cutoff);
vector<AtomPair> allPairsPeriodic(vector<vec3> const& pos, double cutoff, vec3 boxLengths);
vector<AtomPair> allPairsPeriodic_2(vector<vec3> const& pos, double cutoff, vec3 boxLengths);

class Potential {
protected:
    vector<int> colors_to_groups(vector<int> const& colors);
    
public:
    
    unsigned int dimension;
    
    // Build Hamiltonian matrix
    arma::cx_mat build_Hamiltonian(int numAtoms, vector<AtomPair> const& pairs);
    
    
    // Classical pair potential
    double pair_energy(vector<AtomPair> const& pairs);
    
    // Electronic and classical force
    void force(arma::cx_mat const& dE_dH, vector<AtomPair> const& pairs, vector<vec3>& forces, double& pressure);
    //void force_elec(arma::cx_mat const& dE_dH, vector<AtomPair> const& pairs, vector<vec3>& forces, double& pressure);
    
    virtual int numSpins() = 0;
    virtual double rcut() = 0;
    virtual int numOrbitalsPerSite() = 0;
    
    virtual double phi(double r) = 0;
    virtual double dphi_dr(double r) = 0;
    virtual void fill_TB_hoppings(AtomPair pair,
                                  arma::cx_mat& h,
                                  arma::cx_mat& dh_dx,
                                  arma::cx_mat& dh_dy,
                                  arma::cx_mat& dh_dz) = 0;
};

class ToyModelSOribtal : public Potential {
public:
    
    // ss hopping
    static constexpr double h0 = 24;   // should be positive !!!
    static constexpr double alpha = 1.9;
    
    // pair potential
    static constexpr double v0 = 100;
    static constexpr double beta = 3.53;
    
    static constexpr double rmax = 5.0;
    
    ToyModelSOribtal() { dimension = 3; };
    int numSpins() {return 2;};
    double rcut() {return rmax;}
    int numOrbitalsPerSite() {return 1;}
    
    vec3 dt_dr(vec3 delta);
    double phi(double r);
    double dphi_dr(double r);
    void fill_TB_hoppings(AtomPair pair,
                          arma::cx_mat& h,
                          arma::cx_mat& dh_dx,
                          arma::cx_mat& dh_dy,
                          arma::cx_mat& dh_dz);
    
};

class TI_2D_atom : public Potential {
public:

    // ss hopping
    static constexpr double h0 = 24;   // should be positive !!!
    static constexpr double alpha = 1.9;
    
    // pair potential
    static constexpr double v0 = 100;
    static constexpr double beta = 3.53;

    double M = -0.5, lambda = 0.5, t2 = 0.25;
    
    static constexpr double rmax = 5.0;
    
    TI_2D_atom() { dimension = 2; };
    int numSpins() {return 1;};
    double rcut() {return rmax;}
    int numOrbitalsPerSite() {return 2;}
    
    void set_parameters(double _M, double x_lmd, double _t2) {
        M = _M;
        lambda = x_lmd;
        t2 = _t2;
    }
    
    vec3 dt_dr(vec3 delta);
    double phi(double r);
    double dphi_dr(double r);
    void fill_TB_hoppings(AtomPair pair,
                          arma::cx_mat& h,
                          arma::cx_mat& dh_dx,
                          arma::cx_mat& dh_dy,
                          arma::cx_mat& dh_dz);

};

class TD_2D_atom : public Potential 
{
    public:
    // ss hopping
    static constexpr double h0 = 24;   // should be positive !!!
    static constexpr double alpha = 1.9;
    
    // pair potential
    static constexpr double v0 = 100;
    static constexpr double beta = 3.53;

    double M = -0.5, lambda = 0.5, t2 = 0.25;
    
    static constexpr double rmax = 4.0;
    
    TD_2D_atom() { dimension = 2; };
    int numSpins() {return 1;};
    double rcut() {return rmax;}
    int numOrbitalsPerSite() {return 2;}
    
    void set_parameters(double _M, double x_lmd, double _t2) {
        M = _M;
        lambda = x_lmd;
        t2 = _t2;
    }
    
    vec3 dt_dr(vec3 delta);
    double phi(double r);
    double dphi_dr(double r);
    void fill_TB_hoppings(AtomPair pair,
                          arma::cx_mat& h,
                          arma::cx_mat& dh_dx,
                          arma::cx_mat& dh_dy,
                          arma::cx_mat& dh_dz);

};


#endif /* potential_hpp */

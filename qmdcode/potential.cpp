//
//  potential.cpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include "potential.hpp"
#include <armadillo>

SplicePolynomial::SplicePolynomial(double f0, double df0, double x0, double x1, double g1, double dg1) {
    this->x0 = x0;
    a0 = f0;
    a1 = df0;
    double dx = x1 - x0;
    double c1 = a0 + a1*dx - g1;
    double c2 = dx*dx;
    double c3 = dx*dx*dx;
    double d1 = a1 - dg1;
    double d2 = 2*dx;
    double d3 = 3*dx*dx;
    double det = c2*d3 - c3*d2;
    a2 = (c3*d1 - c1*d3) / det;
    a3 = (c1*d2 - c2*d1) / det;
}

double SplicePolynomial::eval(double x) {
    double dx = x - x0;
    return a0 + dx*(a1 + dx*(a2 + dx*a3));
}

double SplicePolynomial::deriv(double x) {
    double dx = x - x0;
    return a1 + dx*(2*a2 + dx*3*a3);
}


std::ostream& operator<< (std::ostream& os, AtomPair const& p) {
    os << "AtomPair { " << p.index1 << ", " << p.index2 << ", " << p.delta << " }";
    return os;
}

vector<AtomPair> allPairs(vector<vec3> const& pos, double cutoff) {
    vector<AtomPair> pairs;
    for (int i = 0; i < pos.size(); i++) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = pos[j]-pos[i];
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {i, j, delta});
        }
    }
    return pairs;
}

/*
double wrapDelta1d(double dx, double boxLength) {
    while (dx >= boxLength/2.0)
        dx -= boxLength;
    while (dx < -boxLength/2.0)
        dx += boxLength;
    return dx;
}
*/

double wrapDelta1d(double dx, double boxLength)
{
    double x_min = -boxLength/2.0;
    double x_max = boxLength/2.0;
    if (dx < x_min)
        dx = x_max - fmod((x_min - dx),(x_max - x_min));
    else
        dx = x_min + fmod((dx - x_min),(x_max - x_min));
    return dx;
}

vec3 wrapDelta(vec3 delta, vec3 boxLengths) {
    delta.x = wrapDelta1d(delta.x, boxLengths.x);
    delta.y = wrapDelta1d(delta.y, boxLengths.y);
    delta.z = wrapDelta1d(delta.z, boxLengths.z);
    return delta;
}

/*double wrapPosition1d(double x, double bdsLo, double bdsHi) {
    double length = bdsHi - bdsLo;
    while (x > bdsHi) {
        x -= length;
    }
    while (x < bdsLo) {
        x += length;
    }
    return x;
}*/

double wrapPosition1d(double x, double bdsLo, double bdsHi) {
    double x_min = bdsLo;
    double x_max = bdsHi;
    if (x < x_min)
        x = x_max - fmod((x_min - x),(x_max - x_min));
    else
        x = x_min + fmod((x - x_min),(x_max - x_min));
    return x;
}

vec3 wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi) {
    p.x = wrapPosition1d(p.x, bdsLo.x, bdsHi.x);
    p.y = wrapPosition1d(p.y, bdsLo.y, bdsHi.y);
    p.z = wrapPosition1d(p.z, bdsLo.z, bdsHi.z);
    return p;
}

vector<AtomPair> allPairsPeriodic(vector<vec3> const& pos, double cutoff, vec3 boxLengths) {
#ifdef WITH_TBB
    tbb::concurrent_vector<AtomPair> pairs;
    auto fn = [&](size_t i) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = wrapDelta(pos[j]-pos[i], boxLengths);
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {(int)i, j, delta});
        }
    };
    tbb::parallel_for(size_t(0), size_t(pos.size()), fn);
    return Vec<AtomPair>(pairs.begin(), pairs.end());
#else
    vector<AtomPair> pairs;
    for (int i = 0; i < pos.size(); i++) {
        for (int j = i+1; j < pos.size(); j++) {
            vec3 delta = wrapDelta(pos[j]-pos[i], boxLengths);
            if (delta.norm2() < cutoff*cutoff)
                pairs.push_back(AtomPair {i, j, delta});
        }
    }
    return pairs;
#endif
}

arma::cx_mat Potential::build_Hamiltonian(int numAtoms, vector<AtomPair> const& pairs) {
    
    int numOrbs = numOrbitalsPerSite();
    arma::cx_mat h(numOrbs, numOrbs);
    arma::cx_mat tmp(numOrbs, numOrbs);
    
    int dim = numAtoms*numOrbs;
    
    arma::cx_mat H(dim, dim);
    H.zeros();
    
    // same site orbital overlap
    for (int i = 0; i < numAtoms; i++) {
        fill_TB_hoppings({i, i, vec3{0,0,0}}, h, tmp, tmp, tmp);
        for (int o1 = 0; o1 < numOrbs; o1++) {
            for (int o2 = 0; o2 < numOrbs; o2++) {
                auto v = h(o1, o2);
                // ensure diagonal is non-zero
                if (o1 == o2 && std::abs(v) == 0.0) {
                    v = std::numeric_limits<double>::epsilon();
                }
                H(i*numOrbs+o1, i*numOrbs+o2) = v;
            }
        }
    }
    
    // neighbor overlap
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            assert(i < j);
            fill_TB_hoppings(p, h, tmp, tmp, tmp);
            for (int o1 = 0; o1 < numOrbs; o1++) {
                for (int o2 = 0; o2 < numOrbs; o2++) {
                    
                    auto v = h(o1, o2);
                    
                    H(i*numOrbs+o1, j*numOrbs+o2) = v;
                    H(j*numOrbs+o2, i*numOrbs+o1) = conj(v);
                }
            }
        }
    }
    
    return H;
}

double Potential::pair_energy(vector<AtomPair> const& pairs) {
    double E_pair = 0.0;
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            E_pair += phi(p.delta.norm());
        }
    }
    return E_pair;
}


void Potential::force(arma::cx_mat const& dE_dH, vector<AtomPair> const& pairs, vector<vec3>& forces, double &virial) {
    std::fill(forces.begin(), forces.end(), vec3{0, 0, 0});
    
    virial = 0;
    
    int nOrbs = numOrbitalsPerSite();
    arma::cx_mat tmp(nOrbs, nOrbs);
    arma::cx_mat dh_dx(nOrbs, nOrbs);
    arma::cx_mat dh_dy(nOrbs, nOrbs);
    arma::cx_mat dh_dz(nOrbs, nOrbs);
    
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            vec3 f_ij {0, 0, 0};
            
            // electronic part
            fill_TB_hoppings(p, tmp, dh_dx, dh_dy, dh_dz);
            for (int o1 = 0; o1 < nOrbs; o1++) {
                for (int o2 = 0; o2 < nOrbs; o2++) {
                    cx_double dE_dH_ij = dE_dH(i*nOrbs+o2, j*nOrbs+o1);
                    cx_double dE_dH_ji = (i == j) ? 0.0 : dE_dH(j*nOrbs+o1, i*nOrbs+o2);
//                    cx_double dE_dH_ji = conj(dE_dH_ij);

// wrong version:
//                    f_ij.x += -numSpins() * real(dh_dx(o2, o1) * dE_dH_ij + dh_dx(o2, o1) * dE_dH_ji);
//                    f_ij.y += -numSpins() * real(dh_dy(o2, o1) * dE_dH_ij + dh_dy(o2, o1) * dE_dH_ji);
//                    f_ij.z += -numSpins() * real(dh_dz(o2, o1) * dE_dH_ij + dh_dz(o2, o1) * dE_dH_ji);

                    f_ij.x += -numSpins() * real(conj(dh_dx(o2, o1)) * dE_dH_ij + dh_dx(o2, o1) * dE_dH_ji);
                    f_ij.y += -numSpins() * real(conj(dh_dy(o2, o1)) * dE_dH_ij + dh_dy(o2, o1) * dE_dH_ji);
                    f_ij.z += -numSpins() * real(conj(dh_dz(o2, o1)) * dE_dH_ij + dh_dz(o2, o1) * dE_dH_ji);

                }
            }
            
            // classical part
            double r = p.delta.norm();
            f_ij += p.delta * (-dphi_dr(r) / r);
            
            forces[j] += f_ij;
            forces[i] -= f_ij;
            virial += f_ij.dot(p.delta);
        }
    }
    
    if(dimension == 2) {
        for(int i=0; i<forces.size(); i++) forces[i].z = 0;
    }
}

/*
void Potential::force_elec(arma::cx_mat const& dE_dH, vector<AtomPair> const& pairs, vector<vec3>& forces, double &virial) {
    std::fill(forces.begin(), forces.end(), vec3{0, 0, 0});
    
    virial = 0;
    
    int nOrbs = numOrbitalsPerSite();
    arma::cx_mat tmp(nOrbs, nOrbs);
    arma::cx_mat dh_dx(nOrbs, nOrbs);
    arma::cx_mat dh_dy(nOrbs, nOrbs);
    arma::cx_mat dh_dz(nOrbs, nOrbs);
    
    for (auto const& p : pairs) {
        if (p.delta.norm2() < rcut()*rcut()) {
            int i = p.index1;
            int j = p.index2;
            vec3 f_ij {0, 0, 0};
            
            // electronic part
            fill_TB_hoppings(p, tmp, dh_dx, dh_dy, dh_dz);
            for (int o1 = 0; o1 < nOrbs; o1++) {
                for (int o2 = 0; o2 < nOrbs; o2++) {
                    cx_double dE_dH_ij = dE_dH(i*nOrbs+o1, j*nOrbs+o2);
                    cx_double dE_dH_ji = (i == j) ? 0.0 : dE_dH(j*nOrbs+o2, i*nOrbs+o1);

                    f_ij.x += -numSpins() * real(dh_dx(o1, o2) * dE_dH_ij + conj(dh_dx(o1, o2)) * dE_dH_ji);
                    f_ij.y += -numSpins() * real(dh_dy(o1, o2) * dE_dH_ij + conj(dh_dy(o1, o2)) * dE_dH_ji);
                    f_ij.z += -numSpins() * real(dh_dz(o1, o2) * dE_dH_ij + conj(dh_dz(o1, o2)) * dE_dH_ji);
                }
            }
            
            forces[j] += f_ij;
            forces[i] -= f_ij;
            virial += f_ij.dot(p.delta);
        }
    }
    
    if(dimension == 2) {
        for(int i=0; i<forces.size(); i++) forces[i].z = 0;
    }

}
*/

//===================================================================================================
// codes for toy-model s-orbital
//===================================================================================================

double ToyModelSOribtal::phi(double r) {
    
    return (r > rmax) ? 0 : v0 * exp(-beta * r);
}

double ToyModelSOribtal::dphi_dr(double r) {
    return (r > rmax) ? 0 : -v0 * beta * exp(-beta * r);
}


vec3 ToyModelSOribtal::dt_dr(vec3 delta) {
    double r = delta.norm();
    
    return (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * delta / r;
}

void ToyModelSOribtal::fill_TB_hoppings(AtomPair pair,
                                        arma::cx_mat& h,
                                        arma::cx_mat& dh_dx,
                                        arma::cx_mat& dh_dy,
                                        arma::cx_mat& dh_dz) {
    h.fill(0);
    dh_dx.fill(0);
    dh_dy.fill(0);
    dh_dz.fill(0);
    if (pair.index1 == pair.index2) {
        h(0, 0) = 0;
        return;
    }
    
    double r = pair.delta.norm();
    
    // ss hopping
    double t     = (r > rmax) ? 0 : -h0 * exp(-alpha * r);
    vec3   dt_dr = (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * pair.delta / r;
    
    h(0, 0) = t;
    dh_dx(0, 0) = dt_dr.x;
    dh_dy(0, 0) = dt_dr.y;
    dh_dz(0, 0) = dt_dr.z;
}

//===================================================================================================
// codes for TI-2D-atom-model
//===================================================================================================

double TI_2D_atom::phi(double r) {
    
    return (r > rmax) ? 0 : v0 * exp(-beta * r);
}

double TI_2D_atom::dphi_dr(double r) {
    return (r > rmax) ? 0 : -v0 * beta * exp(-beta * r);
}

vec3 TI_2D_atom::dt_dr(vec3 delta) {
    double r = delta.norm();
    
    return (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * delta / r;
}

void TI_2D_atom::fill_TB_hoppings(AtomPair pair,
                                        arma::cx_mat& h,
                                        arma::cx_mat& dh_dx,
                                        arma::cx_mat& dh_dy,
                                        arma::cx_mat& dh_dz) {
    h.fill(0);
    dh_dx.fill(0);
    dh_dy.fill(0);
    dh_dz.fill(0);
    
    if (pair.index1 == pair.index2) {
        h(0, 0) = 2. + M;
        h(1, 1) = -(2. + M);
        h(0, 1) = (1. - _I) * lambda;
        h(1, 0) = (1. + _I) * lambda;
        return;
    }
    
    double r = pair.delta.norm();
    auto dlt = pair.delta;
    
    // ss hopping
    double t     = (r > rmax) ? 0 : -h0 * exp(-alpha * r);
    vec3   dt_dr = (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * pair.delta / r;
    
        
    h(0, 0) = 0.5 * t * (-1. + t2);
    h(1, 1) = 0.5 * t * ( 1. + t2);
    
    h(0, 1) = 0.5 * t * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.));
    h(1, 0) = 0.5 * t * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.));
    
    
    dh_dx(0, 0) = 0.5 * dt_dr.x * (-1. + t2);
    dh_dx(1, 1) = 0.5 * dt_dr.x * ( 1. + t2);
    
    dh_dx(0, 1) = 0.5 * dt_dr.x * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-_I/r + _I * dlt.x * (dlt.x - _I * dlt.y) / pow(r, 3)
                             - 2. * lambda * (1. + _I) * dlt.x * pow(dlt.y, 2) / pow(r, 4));
    
    dh_dx(1, 0) = 0.5 * dt_dr.x * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-_I/r + _I * dlt.x * (dlt.x + _I * dlt.y) / pow(r, 3)
                             - 2. * lambda * (1. - _I) * dlt.x * pow(dlt.y, 2) / pow(r, 4));

    
    
    dh_dy(0, 0) = 0.5 * dt_dr.y * (-1. + t2);
    dh_dy(1, 1) = 0.5 * dt_dr.y * ( 1. + t2);
    
    dh_dy(0, 1) = 0.5 * dt_dr.y * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-1/r + _I * dlt.y * (dlt.x - _I * dlt.y) / pow(r, 3)
                             + 2. * lambda * (1. + _I) * (dlt.y / pow(r, 2)) * (1. - pow(dlt.y / r, 2)));
    
    dh_dy(1, 0) = 0.5 * dt_dr.y * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * ( 1/r + _I * dlt.y * (dlt.x + _I * dlt.y) / pow(r, 3)
                             + 2. * lambda * (1. - _I) * (dlt.y / pow(r, 2)) * (1. - pow(dlt.y / r, 2)));

    

    
    // test:
    /*
    double lx = 0.8;
    
    h(0, 0) = 0.5 * t * (-1. + t2);
    h(1, 1) = 0.5 * t * ( 1. + t2);
    
    h(0, 1) = 0.5 * t * (1. + lx * _I);
    h(1, 0) = 0.5 * t * (1. - lx * _I);
    
    
    dh_dx(0, 0) = 0.5 * dt_dr.x * (-1. + t2);
    dh_dx(1, 1) = 0.5 * dt_dr.x * ( 1. + t2);
    
    dh_dx(0, 1) = 0.5 * dt_dr.x * (1. + lx * _I);
    dh_dx(1, 0) = 0.5 * dt_dr.x * (1. - lx * _I);
    
    
    dh_dy(0, 0) = 0.5 * dt_dr.y * (-1. + t2);
    dh_dy(1, 1) = 0.5 * dt_dr.y * ( 1. + t2);
    
    dh_dy(0, 1) = 0.5 * dt_dr.y * (1. + lx * _I);
    dh_dy(1, 0) = 0.5 * dt_dr.y * (1. - lx * _I);
    */
     

}

//===================================================================================================
// codes for TD-2D-atom-model
//===================================================================================================

double TD_2D_atom::phi(double r) {
    
    return (r > rmax) ? 0 : v0 * exp(-beta * r);
}

double TD_2D_atom::dphi_dr(double r) {
    return (r > rmax) ? 0 : -v0 * beta * exp(-beta * r);
}

vec3 TD_2D_atom::dt_dr(vec3 delta) {
    double r = delta.norm();
    
    return (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * delta / r;
}

void TD_2D_atom::fill_TB_hoppings(AtomPair pair,
                                        arma::cx_mat& h,
                                        arma::cx_mat& dh_dx,
                                        arma::cx_mat& dh_dy,
                                        arma::cx_mat& dh_dz) {
    h.fill(0);
    dh_dx.fill(0);
    dh_dy.fill(0);
    dh_dz.fill(0);
    
    if (pair.index1 == pair.index2) {
        h(0, 0) = 2. + M;
        h(1, 1) = -(2. + M);
        h(0, 1) = (1. - _I) * lambda;
        h(1, 0) = (1. + _I) * lambda;
        return;
    }
    
    double r = pair.delta.norm();
    auto dlt = pair.delta;
    
    // ss hopping
    double t     = (r > rmax) ? 0 : -h0 * exp(-alpha * r);
    vec3   dt_dr = (r > rmax) ? vec3(0, 0, 0) : h0 * alpha * exp(-alpha * r) * pair.delta / r;
    
        
    h(0, 0) = 0.5 * t * (-1. + t2);
    h(1, 1) = 0.5 * t * ( 1. + t2);
    
    h(0, 1) = 0.5 * t * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.));
    h(1, 0) = 0.5 * t * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.));
    
    
    dh_dx(0, 0) = 0.5 * dt_dr.x * (-1. + t2);
    dh_dx(1, 1) = 0.5 * dt_dr.x * ( 1. + t2);
    
    dh_dx(0, 1) = 0.5 * dt_dr.x * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-_I/r + _I * dlt.x * (dlt.x - _I * dlt.y) / pow(r, 3)
                             - 2. * lambda * (1. + _I) * dlt.x * pow(dlt.y, 2) / pow(r, 4));
    
    dh_dx(1, 0) = 0.5 * dt_dr.x * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-_I/r + _I * dlt.x * (dlt.x + _I * dlt.y) / pow(r, 3)
                             - 2. * lambda * (1. - _I) * dlt.x * pow(dlt.y, 2) / pow(r, 4));

    
    
    dh_dy(0, 0) = 0.5 * dt_dr.y * (-1. + t2);
    dh_dy(1, 1) = 0.5 * dt_dr.y * ( 1. + t2);
    
    dh_dy(0, 1) = 0.5 * dt_dr.y * (-_I * (dlt.x - _I * dlt.y) / r + lambda * ((1. + _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * (-1/r + _I * dlt.y * (dlt.x - _I * dlt.y) / pow(r, 3)
                             + 2. * lambda * (1. + _I) * (dlt.y / pow(r, 2)) * (1. - pow(dlt.y / r, 2)));
    
    dh_dy(1, 0) = 0.5 * dt_dr.y * (-_I * (dlt.x + _I * dlt.y) / r + lambda * ((1. - _I) * pow(dlt.y / r, 2) - 1.))
                + 0.5 * t * ( 1/r + _I * dlt.y * (dlt.x + _I * dlt.y) / pow(r, 3)
                             + 2. * lambda * (1. - _I) * (dlt.y / pow(r, 2)) * (1. - pow(dlt.y / r, 2)));

    

    
    // test:
    /*
    double lx = 0.8;
    
    h(0, 0) = 0.5 * t * (-1. + t2);
    h(1, 1) = 0.5 * t * ( 1. + t2);
    
    h(0, 1) = 0.5 * t * (1. + lx * _I);
    h(1, 0) = 0.5 * t * (1. - lx * _I);
    
    
    dh_dx(0, 0) = 0.5 * dt_dr.x * (-1. + t2);
    dh_dx(1, 1) = 0.5 * dt_dr.x * ( 1. + t2);
    
    dh_dx(0, 1) = 0.5 * dt_dr.x * (1. + lx * _I);
    dh_dx(1, 0) = 0.5 * dt_dr.x * (1. - lx * _I);
    
    
    dh_dy(0, 0) = 0.5 * dt_dr.y * (-1. + t2);
    dh_dy(1, 1) = 0.5 * dt_dr.y * ( 1. + t2);
    
    dh_dy(0, 1) = 0.5 * dt_dr.y * (1. + lx * _I);
    dh_dy(1, 0) = 0.5 * dt_dr.y * (1. - lx * _I);
    */
     

}
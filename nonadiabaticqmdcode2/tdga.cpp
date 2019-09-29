//
//  model.cpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/2/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include "util.hpp"
#include "model.hpp"
#include "analysis.hpp"

using namespace std;

// {s1, s2} components of pauli matrix vector,
// sigma1     sigma2     sigma3
//  0  1       0 -I       1  0
//  1  0       I  0       0 -1
const Vec3<cx_double> pauli[2][2] {
    {{0, 0, 1}, {1, -_I, 0}},
    {{1, _I, 0}, {0, 0, -1}}
};


void TB_model::init_pp_matrices(void) {
    Mpp = Npp = Upp = arma::mat(3, 3);
    Mpp.zeros();
    Npp.zeros();
    Upp.zeros();
    
    Mpp(0, 1) = Mpp(1, 0) = Mpp(1, 2) = Mpp(2, 1) = sqrt(2.);
    Npp(1, 1) = 0.5;
    Npp(2, 2) = 1;
    Upp(2, 2) = 1;
}


void TB_model::init_random(void) {
    
    //RNG rng = RNG(seed());
    
    std::uniform_real_distribution<double> rd(0, L);
    
    for(int i=0; i<Ns; i++) {
        atom[i].position.x = rd(rng);
        atom[i].position.y = rd(rng);
        atom[i].position.z = rd(rng);
    }
}

void TB_model::compute_hopping_matrix(void) {
    Hopping.zeros();
    for(int i=0; i<Ns; i++) {
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;
//            vec3 delta = (*k).delta;
            double dr = (*k).dr;
            
            assert(dr < r_max);
            Hopping(i, j) = hopping_func(dr);
        }
    }
}

void TB_model::compute_hopping_matrix_full(void) {
    Hopping.zeros();
    for(int i=0; i<Ns; i++) {
        for(int j=0; j<Ns; j++) {
            
            if(j != i) {
                vec3 delta = wrapDelta(atom[j].position - atom[i].position, system_box);
                double dr = delta.norm();
                
                if(dr < r_max)
                    Hopping(i, j) = hopping_func(dr);
                
            }
        }
    }
}

arma::sp_cx_mat TB_model::build_Hamiltonian(void) {  // use atom[i].position to build Hamiltonian (neighbor_list is already updated)
    
    arma::sp_cx_mat H(dim, dim);
    H.zeros();
    
    for(int i=0; i<Ns; i++) {
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;
            double dr = (*k).dr;
            
            H(i, j) = atom[i].R * conj(atom[j].R) * hopping_func(dr);
        }
        
        H(i, i) = atom[i].lambda;
    }
    
    return H;
}

void TB_model::update_atoms(arma::mat & Xm) {
    
    for(int i=0; i<Ns; i++) {
        atom[i].position.x = Xm(0, i);
        atom[i].position.y = Xm(1, i);
        atom[i].position.z = Xm(2, i);
    }
    
    move_atoms_pbc();
    build_cell_atoms_list();
    compute_neighbor_lists();
}

arma::sp_cx_mat TB_model::build_Hamiltonian(arma::mat & Xm) {   // use Xm[ ] to build Hamiltonian
    
    arma::sp_cx_mat H(dim, dim);
    H.zeros();

    for(int i=0; i<Ns; i++) {
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;
                        
            vec3 delta(Xm(0, j) - Xm(0, i), Xm(1, j) - Xm(1, i), Xm(2, j) - Xm(2, i));
            delta = wrapDelta(delta, system_box);

            double dr = delta.norm();
            
            H(i, j) = atom[i].R * conj(atom[j].R) * hopping_func(dr);
        }
        
        H(i, i) = atom[i].lambda;
    }
    
    
    return H;
}


//arma::sp_cx_mat TB_model::build_Hamiltonian(arma::cx_mat & Dm, arma::cx_mat & f) {
//    compute_renormalizations(Dm, f);
//    arma::sp_cx_mat H = build_Hamiltonian();
//    return H;
//}

void TB_model::compute_renormalizations(arma::cx_mat & Dm, arma::cx_mat & f) {
    
    for(int i=0; i<Ns; i++) {
        double n0 = real(Dm(i, i));
        //double n0 = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));
        
        atom[i].R = (conj(f(2, i)) * f(1, i) + conj(f(1, i)) * f(0, i)) / sqrt(n0 * (1. - n0));
        
        //        double mix = 0.5;
        //        atom[i].R = mix * atom[i].R + (1.-mix) * (conj(f(2, i)) * f(1, i) + conj(f(1, i)) * f(0, i)) / sqrt(n0 * (1. - n0));
        
    }
}

void TB_model::compute_Deltas(arma::cx_mat & Dm) {
    for(int i=0; i<Ns; i++) {
        atom[i].Delta = 0;
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;
            atom[i].Delta += Hopping(i, j) * conj(atom[j].R) * Dm(j, i);
        }
    }
}

void TB_model::compute_Deltas(arma::cx_mat & Dm, arma::mat & Xm) {
    for(int i=0; i<Ns; i++) {
        atom[i].Delta = 0;
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            
            int j = (*k).idx;

            vec3 delta(Xm(0, j) - Xm(0, i), Xm(1, j) - Xm(1, i), Xm(2, j) - Xm(2, i));
            delta = wrapDelta(delta, system_box);

            double dr = delta.norm();
            
            atom[i].Delta += hopping_func(dr) * conj(atom[j].R) * Dm(j, i);
        }
    }
}


void TB_model::compute_fermi_level(arma::vec & eigE) {
    
    double x1 = eigE(0);
    double x2 = eigE(eigE.size()-1);
    
    int max_bisection = 500;
    double eps_bisection = 1.e-12;
    
    int iter = 0;
    while(iter < max_bisection || fabs(x2 - x1) > eps_bisection) {
        
        double xm = 0.5*(x1 + x2);
        double density = 0;
        for(unsigned int i=0; i<eigE.size(); i++) {
            density += fermi_density(eigE(i), kT, xm);
        }
        density /= ((double) dim);
        
        if(density <= filling) x1 = xm;
        else x2 = xm;
        
        iter++;
    }
    
    mu = 0.5*(x1 + x2);
}

arma::cx_mat TB_model::compute_density_matrix(arma::sp_cx_mat & H) {
    arma::cx_mat Hm(H);
    return compute_density_matrix(Hm);
}

arma::cx_mat TB_model::compute_density_matrix(arma::cx_mat & H) {
    
    arma::cx_mat Hm(H);
    arma::cx_mat D(dim, dim);
    
    arma::vec eigval(dim);
    arma::cx_mat eigvec;
    
    arma::eig_sym(eigval, eigvec, Hm);
    
    compute_fermi_level(eigval);
    
    arma::vec fd_factor(dim);
    for(int i=0; i<dim; i++) fd_factor(i) = fermi_density(eigval(i), kT, mu);
    
    D.zeros();
    for(int a=0; a<dim; a++)
        for(int b=a; b<dim; b++) {
            
            cx_double sum = 0;
            for(int m=0; m<dim; m++) {
                sum += fd_factor(m) * conj(eigvec(a, m)) * eigvec(b, m);
            }
            D(a, b) = sum;
            if(a != b) D(b, a) = conj(sum);
        }
    
    return D;
}

arma::cx_mat TB_model::compute_tV_density_matrix(arma::sp_cx_mat & H) {
    
    arma::cx_mat Hm(H);
    arma::cx_mat D(dim, dim);
    
    arma::vec eigval(dim);
    arma::cx_mat eigvec;
    
    arma::eig_sym(eigval, eigvec, Hm);
    
    compute_fermi_level(eigval);
    
    arma::vec fd_factor(dim);
    for(int i=0; i<dim; i++) fd_factor(i) = fermi_density(eigval(i), kT, mu);
    
    D.zeros();
    for(int i=0; i<Ns; i++) {
        
        cx_double sum = 0;
        for(int m=0; m<dim; m++) {
            sum += fd_factor(m) * conj(eigvec(i, m)) * eigvec(i, m);
        }
        D(i, i) = sum;
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            int j = (*k).idx;
            
            if(j > i) {
                
                cx_double sum2 = 0;
                for(int m=0; m<dim; m++) {
                    sum2 += fd_factor(m) * conj(eigvec(j, m)) * eigvec(i, m);
                }
                D(i, j) = sum2;
                D(j, i) = conj(sum2);
                
            }
        }
        
    }
    return D;
}

arma::cx_mat TB_model::init_uncorrelated_Phi(arma::cx_mat & Dm) {
    
    arma::cx_mat f(3, dim);
    for(int i=0; i<Ns; i++) {
        double n0 = real(Dm(i, i));
        f(0, i) = 1. - n0;
        f(1, i) = sqrt(n0 * (1. - n0));
        f(2, i) = n0;
        
//        f(0, i) = 0.5;
//        f(1, i) = 1./sqrt(2.);
//        f(2, i) = 0.5;
        
    }
    
    return f;
}

arma::cx_mat TB_model::solve_Phi(arma::cx_mat & Dm) {
    
    //compute_renormalizations(Dm, Phi);
    
    compute_Deltas(Dm);
    
    arma::cx_mat f(3, dim);
    
    for(int i=0; i<Ns; i++) {
        arma::cx_mat Hpp(3, 3);
        double n0 = real(Dm(i, i));
        //double n0 = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));
        
        auto tmpx = atom[i].R * atom[i].Delta;
        double Lmd = 2. * real(tmpx) * (n0 - 0.5) / (n0 * (1. - n0));
        
        Hpp = (atom[i].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp + (2. * Lmd - atom[i].lambda) * Npp;
        
        arma::vec E(3);
        arma::cx_mat v(3, 3);
        
        arma::eig_sym(E, v, Hpp);
        
        int idx_min = 0;
        
        double e_min = E(0);
        for(int k=1; k<3; k++) {
            if(E(k) < e_min) {
                e_min = E(k);
                idx_min = k;
            }
        }
        
        f(0, i) = v(0, idx_min);
        f(1, i) = v(1, idx_min) / sqrt(2.);
        f(2, i) = v(2, idx_min);
        
    }
    
    return f;
}

double TB_model::adjust_lambda(arma::cx_mat & Dm, arma::cx_mat & f) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        double n0 = real(Dm(i, i));
        double npp = real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i));
        
        //        int sgn = (n0 > npp) ? +1 : -1;
        //        atom[i].lambda += sgn * 0.1 * pow(fabs(n0 - npp), 0.3); // 0.01 * sgn * pow(fabs(n0 - npp), 0.8);
        
        atom[i].lambda += 0.5 * (n0 - npp);
        
        sum += fabs(n0 - npp);
    }
    
    return sum/((double) Ns);
}

void TB_model::test(void) {
    
    /*
    int ix = 2;
    
    arma::cx_mat Hpp(3, 3);
    double n0 = real(Density_Mat(ix, ix));
    
    Hpp = (2. * atom[ix].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp - atom[ix].lambda * Npp;
    
    arma::vec E(3);
    arma::cx_mat v(3, 3);
    
    arma::eig_sym(E, v, Hpp);
    
    cout << Hpp << endl;
    
    cout << E << endl;
    cout << v << endl;
    */
    
    ofstream fs("tst.dat");
    for(int i=0; i<Ns; i++) {
        int ic = atom[i].cell_idx;
        fs << i << "\t (" << cell[ic].x << ", " << cell[ic].y << ", " << cell[ic].z << ")" << endl;
        fs << "...................." << endl;
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
         
            int j = (*k).idx;
            ic = atom[j].cell_idx;
            
            fs << j << "\t (" << cell[ic].x << ", " << cell[ic].y << ", " << cell[ic].z << ")" << endl;

        }
        fs << endl;
        fs << "------------------------------------" << endl;
    }
    fs.close();
}

void TB_model::compute_GA_solution(void) {
    
    arma::cx_mat Dm(dim, dim);
    arma::cx_mat f(3, dim);
    arma::cx_mat f_prev(3, dim);
    
    kT = GA_kT_start;
    
    double error = 100;
    double err_f = 100;
    int max_iter = 100000;
    int iter = 0;
    
    reset_GA_parameters();
    
    Hamiltonian = build_Hamiltonian();
    
    Dm = compute_tV_density_matrix(Hamiltonian);
    f = init_uncorrelated_Phi(Dm);
    f_prev = f;
    
    double tag_zeroR = 10;
    
    ofstream fx("b.dat");
    
    while(((iter < max_iter && error > 1.e-5) || kT > GA_kT_stop) && tag_zeroR > 1.e-6) {
        
        cout << "iter-" << iter << '\t';
        
        compute_renormalizations(Dm, f);
        Hamiltonian = build_Hamiltonian();
        
        Dm = compute_tV_density_matrix(Hamiltonian);
        compute_Deltas(Dm);
        f = solve_Phi(Dm);
        
        auto df = f - f_prev;
        err_f = arma::norm(df)/(3.*dim);
        f_prev = f;
        
        tag_zeroR = 0;
        for(int i=0; i<Ns; i++) tag_zeroR += real(f(2, i) * conj(f(2, i)));
        tag_zeroR /= ((double) Ns);
        
        error = adjust_lambda(Dm, f);
        
        cout << "err = " << error << '\t' << err_f << '\t' << tag_zeroR << endl;
        
        iter++;
        
        fx << error << '\t' << kT << endl;
        
        if(iter % GA_n_annealing == 0 && kT > GA_kT_stop) {
            kT *= GA_r_annealing;
            cout << "kT = " << kT << endl;
            cout << "err = " << error << '\t' << err_f << '\t' << tag_zeroR << endl;
            
        }
        
        
        if(iter % 20 == 0) {
            ofstream fs("t.dat");
            fs.precision(15);
            for(int i=0; i<Ns; i++) {
                fs << i << '\t' << real(atom[i].R) << '\t' << real(f(0, i)) << '\t' << real(f(1, i)) << '\t' << real(f(2, i)) << '\t' << real(Dm(i, i)) << '\t' << real(conj(f(1, i)) * f(1, i) + conj(f(2, i)) * f(2, i)) << '\t' << atom[i].lambda << endl;
            }
            fs.close();
            
        }
        
        //        ofstream fxx("m.dat", ios::app);
        //        int ix = 5;
        //        arma::cx_mat Hpp(3, 3);
        //        double n0 = real(Dm(ix, ix));
        //        cout << "n0 = " << n0 << endl;
        //        compute_Deltas(Dm);
        //        double e_tnn = 0;
        //        cx_double _delta = 0;
        //        for(int k=0; k<N_nn1; k++) {
        //            int j = atom[ix].nn1[k]->idx;
        //            e_tnn += t1 * real(Dm(j, ix));
        //            _delta += t1 * conj(atom[j].R) * Dm(j, ix);
        //        }
        //        cout << "e_k = " << e_tnn << "\t R = " << real(atom[ix].R) << "\t Dlt = " << _delta.real() << "," << _delta.imag() << endl;
        //
        //        Hpp = (atom[ix].Delta / sqrt(n0 * (1. - n0))) * Mpp + U * Upp - atom[ix].lambda * Npp;
        //        cout << Hpp << endl;
        //        fxx.close();
        
    }
    
    fx.close();
    
    
    cout << "computing full density matrix ..." << endl;
    Density_Mat = compute_density_matrix(Hamiltonian);
    Phi = f;
    
}

void TB_model::analyze_data(void) {
    
    double sum1_nd = 0, sum2_nd = 0;
    double sum1_np = 0, sum2_np = 0;
    double sum1_R = 0, sum2_R = 0;
    double sum1_d = 0, sum2_d = 0;
    
    for(int i=0; i<Ns; i++) {
        
        double nr = real(Density_Mat(i, i));
        sum1_nd += nr;
        sum2_nd += pow(nr, 2);
        
        double npp = real(conj(Phi(1, i)) * Phi(1, i) + conj(Phi(2, i)) * Phi(2, i));
        sum1_np += npp;
        sum2_np += pow(npp, 2);
        
        double db = real(conj(Phi(2, i)) * Phi(2, i));
        sum1_d += db;
        sum2_d += pow(db, 2);
        
        double rq = real(atom[i].R * conj(atom[i].R));

        sum1_R += rq;
        sum2_R += pow(rq, 2);
        
    }
    
    sum1_nd /= ((double) Ns);
    sum2_nd /= ((double) Ns);
    sum1_np /= ((double) Ns);
    sum2_np /= ((double) Ns);
    sum1_d /= ((double) Ns);
    sum2_d /= ((double) Ns);
    sum1_R /= ((double) Ns);
    sum2_R /= ((double) Ns);
    
    avg_nd = sum1_nd;
    sgm_nd = sqrt(sum2_nd - sum1_nd * sum1_nd);
    
    avg_np = sum1_np;
    sgm_np = sqrt(sum2_np - sum1_np * sum1_np);
    
    avg_d = sum1_d;
    sgm_d = sqrt(sum2_d - sum1_d * sum1_d);
    
    avg_R = sum1_R;
    sgm_R = sqrt(sum2_R - sum1_R * sum1_R);
    
}

void TB_model::compute_electronic_energy(void) {
    
    cx_double sum = 0;
    cx_double sum2 = 0;
    
    for(int i=0; i<Ns; i++) {
        cx_double tmp;
        
        tmp = Hamiltonian(i, i);
        //        sum += (tmp - atom[i].lambda) * Density_Mat(i, i);
        sum += (tmp) * Density_Mat(i, i);
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
            int j = (*k).idx;
            tmp = Hamiltonian(i, j);
            sum += tmp * Density_Mat(j, i);
        }
        
        sum2 += Phi(2, i) * conj(Phi(2, i));
    }
    
    
    E_e = (2. * sum.real() + U * sum2.real()) / ((double) Ns);
}

void TB_model::compute_pair_potential(void) {
    
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {
//            int j = (*k).idx;
            double dr = (*k).dr;
            
            sum += phi(dr);
        }
    }
    
    E_p = 0.5 * sum / ((double) Ns);
}

void TB_model::compute_kinetic_energy(void) {
    double sum = 0;
    for(int i=0; i<Ns; i++) {
        vec3 P(momentum(0, i), momentum(1, i), momentum(2, i));
        sum += P.norm2();
    }
    
    E_k = 0.5 * sum / mass / ((double) Ns);
}

arma::cx_mat TB_model::compute_dF(arma::cx_mat & F, arma::cx_mat &Dm) {
    
    arma::cx_mat dF(3, dim);
    
    for(int i=0; i<Ns; i++) {
        
        double n0 = real(Dm(i, i));
        //double n0 = real(conj(F(1, i)) * F(1, i) + conj(F(2, i)) * F(2, i));
        
        double qr = sqrt(n0 * (1. - n0));
        
        auto tmpx = atom[i].R * atom[i].Delta;
        double Lmd = 2. * real(tmpx) * (n0 - 0.5) / (n0 * (1. - n0));
        
        
        dF(0, i) = -_I * 2. * F(1, i) * conj(atom[i].Delta) / qr;
        
        dF(1, i) = -_I * (F(2, i) * conj(atom[i].Delta) + F(0, i) * atom[i].Delta) / qr - _I * Lmd * F(1, i) + 0.5 * _I * atom[i].lambda * F(1, i);
        
        dF(2, i) = -_I * U * F(2, i) - _I * 2. * F(1, i) * atom[i].Delta / qr - 2. * _I * Lmd * F(2, i) + _I * atom[i].lambda * F(2, i);
        
        
        //        for(int k=0; k<N_nn1; k++) {
        //            int j = atom[i].nn1[k]->idx;
        //
        //            cout << Dm(i, j) << endl;
        //        }
        //
        //        cout << dF(0, i) << '\t' << dF(1, i) << '\t' << dF(2, i) << '\t' << U << '\t' << Lmd << endl;
        //        exit(1);
        
    }
    
    return dF;
}

arma::mat TB_model::compute_dX(arma::mat & Pm) {
    arma::mat dX(3, dim);
    for(int i=0; i<Ns; i++) {
        
        for(int k=0; k<3; k++) dX(k, i) = Pm(k, i) / mass;        
    }
    return dX;
}


arma::mat TB_model::compute_dP(arma::mat & Xm, arma::cx_mat & Dm) {
    arma::mat dP(3, dim);
    
    dP.zeros();
    
    for(int i=0; i<Ns; i++) {
        
        for(auto k=atom[i].neighbor_list.begin(); k != atom[i].neighbor_list.end(); ++k) {

            int j = (*k).idx;
            double dlx[3];
            for(int m=0; m<3; m++) dlx[m] = Xm(m, j) - Xm(m, i);
            vec3 delta(dlx[0], dlx[1], dlx[2]);
            
            delta = wrapDelta(delta, system_box);
            
            vec3 fi = dphi_dr(delta) + 4. * real(atom[i].R * conj(atom[j].R) * Dm(j, i)) * dt_dr(delta);
            
            dP(0, i) += fi.x;
            dP(1, i) += fi.y;
            dP(2, i) += fi.z;
            
        }
    }
    
    return dP;
}

void TB_model::integrate_EOM_RK4(double dt) {
    
    arma::cx_mat H(Hamiltonian);
    arma::cx_mat D(Density_Mat);
    arma::cx_mat D2, KD, KD_sum;
    
    arma::cx_mat F(Phi);
    arma::cx_mat F2, KF, KF_sum;
    
    arma::mat X(position);
    arma::mat P(momentum);
    arma::mat X2, KX, KX_sum;
    arma::mat P2, KP, KP_sum;


    // ------- RK4 step-1: ----------------
    
    KD = -_I * dt * ( H * D - D * H );
    KF = dt * compute_dF(F, D);
    KX = dt * compute_dX(P);
    KP = dt * compute_dP(X, D);
    
    D2 = D + 0.5 * KD;
    KD_sum = KD / 6.;
    
    F2 = F + 0.5 * KF;
    KF_sum = KF / 6.;
    
    X2 = X + 0.5 * KX;
    KX_sum = KX / 6.;
    
    P2 = P + 0.5 * KP;
    KP_sum = KP / 6.;

    // ------- RK4 step-2: ----------------
    
    compute_renormalizations(D2, F2);
    update_atoms(X2);
    H = build_Hamiltonian(X2);
    KD = -_I * dt * ( H * D2 - D2 * H );
    
    compute_Deltas(D2, X2);
    KF = dt * compute_dF(F2, D2);
    KX = dt * compute_dX(P2);
    KP = dt * compute_dP(X2, D2);
    
    D2 = D + 0.5 * KD;
    KD_sum += KD / 3.;
    
    F2 = F + 0.5 * KF;
    KF_sum += KF / 3.;
    
    X2 = X + 0.5 * KX;
    KX_sum += KX / 3.;
    
    P2 = P + 0.5 * KP;
    KP_sum += KP / 3.;
    
    // ------- RK4 step-3: ----------------
    
    compute_renormalizations(D2, F2);
    update_atoms(X2);
    H = build_Hamiltonian(X2);
    KD = -_I * dt * ( H * D2 - D2 * H );
    
    compute_Deltas(D2, X2);
    KF = dt * compute_dF(F2, D2);
    KX = dt * compute_dX(P2);
    KP = dt * compute_dP(X2, D2);
    
    D2 = D + KD;
    KD_sum += KD / 3.;
    
    F2 = F + KF;
    KF_sum += KF / 3.;
    
    X2 = X + KX;
    KX_sum += KX / 3.;
    
    P2 = P + KP;
    KP_sum += KP / 3.;

    // ------- RK4 step-4: ----------------
    
    compute_renormalizations(D2, F2);
    update_atoms(X2);
    H = build_Hamiltonian(X2);
    KD = -_I * dt * ( H * D2 - D2 * H );
    KD_sum += KD / 6.;
    
    compute_Deltas(D2, X2);
    KF = dt * compute_dF(F2, D2);
    KF_sum += KF / 6.;
    
    KX = dt * compute_dX(P2);
    KX_sum += KX / 6.;
    
    KP = dt * compute_dP(X2, D2);
    KP_sum += KP / 6.;
    
    // ------- RK4: sum all steps: ------------
    
    Density_Mat = D + KD_sum;
    Phi = F + KF_sum;
    
    position = X + KX_sum;
    momentum = P + KP_sum;
    
    // compute the system Hamiltonian, R, Delta:
    
    compute_renormalizations(Density_Mat, Phi);
    update_atoms(position);
    Hamiltonian = build_Hamiltonian(position);
    
    compute_Deltas(Density_Mat, position);


}

void TB_model::init_config(void) {
    for(int i=0; i<Ns; i++) {
        
        //for(int m=0; m<3; m++) momentum(m, i) = 0;
        
        momentum(0, i) = mass * atom[i].velocity.x;
        momentum(1, i) = mass * atom[i].velocity.y;
        momentum(2, i) = mass * atom[i].velocity.z;
        
        position(0, i) = atom[i].position.x;
        position(1, i) = atom[i].position.y;
        position(2, i) = atom[i].position.z;
    }
    
}

void TB_model::compute_cm_velocity(void) {
    
    V_cm.x = V_cm.y = V_cm.z = 0;
    
    for(int i=0; i<Ns; i++) {
        V_cm.x += momentum(0, i);
        V_cm.y += momentum(1, i);
        V_cm.z += momentum(2, i);
    }
    
    V_cm = V_cm / (mass * Ns);
    
}

void TB_model::adjust_cm_velocity(void) {
    compute_cm_velocity();
    
    for(int i=0; i<Ns; i++) {
        momentum(0, i) -= mass * V_cm.x;
        momentum(1, i) -= mass * V_cm.y;
        momentum(2, i) -= mass * V_cm.z;
    }
}

void TB_model::simulate_quench(int max_steps, double dt, double U_i, double U_f) {
    
    
    //-----------------------------------------------------
    // initialize the elastic sub-system:
    
    
    //-----------------------------------------------------
    
    /*
    init_random();
    build_cell_atoms_list();
    compute_neighbor_lists();
    compute_hopping_matrix();
    */
        
    int nstep_annealing = 10000;
    double dt_annealing = 0.05;
    gamma = 0.1;
    QMD_simulation(nstep_annealing, dt_annealing);

    
    move_atoms_pbc();
    build_cell_atoms_list();
    compute_neighbor_lists();
    compute_hopping_matrix();
    
    
    U = U_i;
    compute_GA_solution();
    
    U = U_f;

    // transfer QMD data to GvNMD simulation
    init_config();
    adjust_cm_velocity();
    
    int n_save = 100;
    
    ofstream fs("a.dat");

    // ===============================================

    MD_Data accu_data;
    MD_Data mddata;
    double r_min = 0.0002;
    double r_max = 0.5 * L;
    double dr = 0.02;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    double D_max = 0.3;
    int num_hs_D = 250;
    double dlt_D = D_max/((double) num_hs_D);
    mddata.init_hist_D(num_hs_D, dlt_D);
    
    double R_max = 1.0;
    int num_hs_R = 250;
    double dlt_R = R_max/((double) num_hs_R);
    mddata.init_hist_R(num_hs_R, dlt_R);
    
    
    int n_measure_steps = 2000;
    int N_measurement = 20;
    int n_separation = 100;
    
    accu_data.init_time_series(n_measure_steps, dt);

    Vec<Measurement> measurement(N_measurement);
    for(int k=0; k<N_measurement; k++) measurement[k] = Measurement(dt, n_measure_steps, Ns);


    // ===============================================

    int nx = 0;
    
    for(int i=0; i<max_steps; i++) {
        
        if(i % 100 == 0)
            cout << "i = " << i << endl;
        
        integrate_EOM_RK4(dt);
        
        double n_tot = real( arma::trace(Density_Mat) ) / ((double) dim);
        analyze_data();
        compute_pair_potential();
        compute_kinetic_energy();
        compute_electronic_energy();
        
        compute_cm_velocity();
        
        fs << i * dt << '\t' << avg_R << '\t' << sgm_R << '\t' << avg_d << '\t' << sgm_d << '\t';
        fs << avg_nd << '\t' << sgm_nd << '\t' << avg_np << '\t' << sgm_np << '\t';
        fs << E_e << '\t' << E_p << '\t' << E_k << '\t';
        fs << V_cm.x << '\t' << V_cm.y << '\t' << V_cm.z << '\t' << n_tot << endl;
        
        if(i % n_save == 0) {
         
            auto pairs = find_all_pairs(1.5 * r_max);
            mddata.compute_g_r(pairs, Ns, L*L*L, r_min, r_max);
            mddata.compute_hist_double_occupancy(*this, D_max);
            mddata.compute_hist_renormalization(*this, R_max);
            mddata.basic_measurements(*this);

//            mddata.print_g_r("gr" + std::to_string(i) + ".dat");
//            mddata.print_hist_D("hd"  + std::to_string(i) + ".dat");
//            mddata.print_hist_R("hr" + std::to_string(i) + ".dat");
            
            if(i > 1000) {
                average_MD_Data(accu_data, mddata, nx);
                nx++;
            }
            accu_data.print_g_r("gr.dat");
            accu_data.print_hist_D("hd.dat");
            accu_data.print_hist_R("hr.dat");

            accu_data.print_data("s.dat", U_f);
        }
        
        for(int k=0; k<N_measurement; k++) {
            
            if(i == (k * n_separation)) measurement[k].start(*this); // the one-time activation
            
            if(measurement[k].active == 1) {
                
                measurement[k].record_data(*this);
                
            } else if(measurement[k].enabled == 1) {
                
                measurement[k].compute_diffusion_coef();
                
                accu_data.average_ms_data(measurement[k]);
                
                accu_data.print_auto_correlation("d.dat");
                
                
                // restart the measurement !!
                measurement[k].start(*this);
                
            }
            
            //cout << measurement[k].active << '\t';
        }

    }
    fs.close();
    
}











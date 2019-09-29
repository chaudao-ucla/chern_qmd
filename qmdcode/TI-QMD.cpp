//
//  main.cpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <iostream>

#include "util.hpp"
#include "potential.hpp"
#include "tbmd.hpp"
#include "analysis.hpp"
#include "test.hpp"


void QMD_simulation(int N_atoms, double density, double kT, int N_steps, double dt) 
{
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./2.);
    cout << "l_box = " << l_box << endl;

    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT, gamma);
    atoms.set_BD_type(1);

    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT = " << kT << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;

    
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0; // 0.5;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;

    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    
    //----- initialize the system ------
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    //atoms.init_from_file(lt, lz, bd_lo, bd_hi);

    
    atoms.init_potential();
    
    MD_Data mddata;
    MD_Data accu_data;
    int n_save = 30;
    int n_reset = 10000;
    int n_bb = 10;
    int nx = 0;
    
    double r_min = 0.0002;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    
    
    ofstream fs("r.dat");
    ofstream fs2("b.dat");
    
    for(int r=0; r < N_steps; r++) {
        
        cout << "r = " << r << endl;
        
        atoms.step_NVT(dt);

        //atoms.test_1();
        
        if(r % n_reset == 0) {
            fs.close();
            fs.open("r.dat");
        }
        
        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();
        double e_tot = ee + ep + ek;
        
        double e_gap = atoms.compute_Eg();
        
        fs << r * dt << '\t' << e_tot<< '\t' << ee << '\t' << ep << '\t' << ek << '\t' << e_gap << '\t';
        fs << endl;
        
        if(r % n_bb == 0) {
            double bott_idx = atoms.compute_Bott_index();
            double sigma_xy = atoms.compute_sigma_xy();
            fs2 << bott_idx << '\t' << sigma_xy << endl;
        }
        
        if(r % n_save == 0) {
            
            mddata.basic_measurements(atoms);
            mddata.compute_g_r(atoms.pairs, atoms.numAtoms, atoms.volume, r_min, r_max, atoms.potential->dimension);
            
            average_MD_Data(accu_data, mddata, nx);
            nx++;
            
            accu_data.print_g_r("gr.dat");
            
            accu_data.print_data("a.dat", density);
            
            atoms.save_configuration("config1.dat");
            
//            if(r % 50 == 0) {
//                ofstream fxx("pos.dat");
//                fxx.precision(10);
//                for(int i=0; i<atoms.numAtoms; i++) {
//                    fxx << atoms.position[i].x << '\t' << atoms.position[i].y << '\t' << atoms.position[i].z << '\t';
//                    fxx << atoms.velocity[i].x << '\t' << atoms.velocity[i].y << '\t' << atoms.velocity[i].z << '\t';
//                    fxx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
//                }
//                fxx.close();
//            }
        }
    }
    
    fs.close();
    fs2.close();
    
    //    atoms.compute_forces_elec();
    //    ofstream fx("ft.dat");
    //    for(int i=0; i<N_atoms; i++) {
    //        fx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
    //    }
    //    fx.close();
    
    
}

void QMD_annealing_simulation(int N_atoms, double density, double kT_init, 
        double r_kT, int N_steps, double dt) 
{
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./2.);
    cout << "l_box = " << l_box << endl;
    
    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT_init, gamma);
    atoms.set_BD_type(1);
    
    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT_init = " << kT_init << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;
    
    
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0; // 0.5;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    atoms.init_potential();
    
    MD_Data mddata;
    MD_Data accu_data;
    int n_save = 30;
    int n_reset = 10000;
    int n_bb = 10;
    int nx = 0;
    
    double r_min = 0.0002;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    
    
    ofstream fs("r.dat");
    ofstream fs2("b.dat");
    
    double time = 0;
    
    for(int r=0; r < N_steps; r++) {
        
        cout << "r = " << r << endl;
    
        atoms.kT = kT_init * exp(-r_kT * time);
    
        atoms.step_NVT(dt);
        
        if(r % n_reset == 0) {
            fs.close();
            fs.open("r.dat");
        }
        
        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();
        
        double e_gap = atoms.compute_Eg();
        
        fs << r * dt << '\t' << ee << '\t' << ep << '\t' << ek << '\t' << e_gap << '\t';
        fs << atoms.kT << '\t';
        fs << endl;
        
        if(r % n_bb == 0) {
            double bott_idx = atoms.compute_Bott_index();
            fs2 << bott_idx << endl;
        }
        
        if(r % n_save == 0) {
            
            mddata.basic_measurements(atoms);
            mddata.compute_g_r(atoms.pairs, atoms.numAtoms, atoms.volume, r_min, r_max, atoms.potential->dimension);
            
            average_MD_Data(accu_data, mddata, nx);
            nx++;
            
            accu_data.print_g_r("gr.dat");
            
            accu_data.print_data("a.dat", density);
            
            atoms.save_configuration("pos.dat");
            
            //            if(r % 50 == 0) {
            //                ofstream fxx("pos.dat");
            //                fxx.precision(10);
            //                for(int i=0; i<atoms.numAtoms; i++) {
            //                    fxx << atoms.position[i].x << '\t' << atoms.position[i].y << '\t' << atoms.position[i].z << '\t';
            //                    fxx << atoms.velocity[i].x << '\t' << atoms.velocity[i].y << '\t' << atoms.velocity[i].z << '\t';
            //                    fxx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
            //                }
            //                fxx.close();
            //            }
        }
        
        time += dt;

    }
    
    fs.close();
    fs2.close();
    
    //    atoms.compute_forces_elec();
    //    ofstream fx("ft.dat");
    //    for(int i=0; i<N_atoms; i++) {
    //        fx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
    //    }
    //    fx.close();
    
    
}

void QMD_simulation_open_BC(int N_atoms, double density, double kT, int N_steps, double dt) 
{
    
    
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./2.);
    cout << "l_box = " << l_box << endl;
    
    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT, gamma);
    atoms.set_BD_type(0);  // <---- IMPORTANT : open BC
    
    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT = " << kT << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;
    
    
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0; // 0.5;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    
    //----- initialize the system ------
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    //atoms.init_from_file(lt, lz, bd_lo, bd_hi);
    
    
    atoms.init_potential();
    
    MD_Data mddata;
    MD_Data accu_data;
    int n_save = 10;
    int n_reset = 10000;
    int n_bb = 10;
    int nx = 0;
    
    double r_min = 0.0002;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    
    
    ofstream fs("r.dat");
    ofstream fs2("b.dat");
    
    for(int r=0; r < N_steps; r++) {
        
        cout << "r = " << r << endl;
        
        atoms.step_NVT(dt);
        
        if(r % n_reset == 0) {
            fs.close();
            fs.open("r.dat");
        }
        
        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();
        
        double e_gap = atoms.compute_Eg();
        
        fs << r * dt << '\t' << ee << '\t' << ep << '\t' << ek << '\t' << e_gap << '\t';
        fs << endl;
        
        /*if(r % n_bb == 0) {
            double bott_idx = atoms.compute_Bott_index();
            fs2 << bott_idx << endl;
        }*/
        
        if(r % n_save == 0) {
            
            mddata.basic_measurements(atoms);
            mddata.compute_g_r(atoms.pairs, atoms.numAtoms, atoms.volume, r_min, r_max, atoms.potential->dimension);
            
            average_MD_Data(accu_data, mddata, nx);
            nx++;
            
            accu_data.print_g_r("gr.dat");
            
            accu_data.print_data("a.dat", density);
            
            atoms.save_configuration("pos.dat");
            
            //            if(r % 50 == 0) {
            //                ofstream fxx("pos.dat");
            //                fxx.precision(10);
            //                for(int i=0; i<atoms.numAtoms; i++) {
            //                    fxx << atoms.position[i].x << '\t' << atoms.position[i].y << '\t' << atoms.position[i].z << '\t';
            //                    fxx << atoms.velocity[i].x << '\t' << atoms.velocity[i].y << '\t' << atoms.velocity[i].z << '\t';
            //                    fxx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
            //                }
            //                fxx.close();
            //            }
            
            int Nq = 20;
            for(int z=0; z<Nq; z++) {
                atoms.plot_eigen_mode("mode" + std::to_string(z) + ".dat", atoms.num_filling + z + 1);
            }
            
            atoms.print_eig_E("eig.dat");
        }
    }
    
    fs.close();
    fs2.close();
    
    //    atoms.compute_forces_elec();
    //    ofstream fx("ft.dat");
    //    for(int i=0; i<N_atoms; i++) {
    //        fx << atoms.force[i].x << '\t' << atoms.force[i].y << '\t' << atoms.force[i].z << endl;
    //    }
    //    fx.close();
    
    
}

void QMD_quench_simulation(int N_atoms, double density, double kT_init,
        double r_kT, int N_steps, double dt, int nstep_annealing = 10000,
        double dt_annealing = 0.05)
{
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./3.);
    cout << "l_box = " << l_box << endl;
    
    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT_init, gamma);
    atoms.set_BD_type(1);
    
    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT_init = " << kT_init << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;
    
    
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0; // 0.5;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    atoms.init_potential();
    
    for(int r=0; r < nstep_annealing; r++)
        atoms.step_NVT(dt_annealing);
    
    for(int i=0; i<N_steps; i++)
    {
        atoms.evolve_NVT_RK4(dt);
    }
}


void QMD_simulation_noneq(int N_atoms, double density, double kT, int N_steps, 
        double dt)
{
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./2.);
    cout << "l_box = " << l_box << endl;


    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT, gamma);

    //atoms.set_BD_type(0);
    atoms.set_BD_type(1);

    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT_init = " << kT << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;
    
    
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0.0;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;

    //atoms.init_from_file(lt, lz, bd_lo, bd_hi);
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    atoms.init_potential();

    MD_Data mddata;
    MD_Data accu_data;
    int n_save = 30;
    int n_reset = 10000;
    int n_bb = 10;
    int nx = 0;
    
    double r_min = 0.0002;
    double r_max = 0.25 * (lt + lz);
    double dr = 0.01;
    int num_r = (int) ((r_max - r_min)/dr);
    mddata.init_g_r(num_r, dr);
    
    ofstream fs("r.dat");
    ofstream fs2("b.dat");
    double squaremean = 0;
    double meansquare = 0; 

    for(int r=0; r < N_steps; r++) {
        
        cout << "r = " << r << endl;
    
        atoms.evolve_NVT_RK4(dt);
        
        if(r % n_reset == 0) {
            fs.close();
            fs.open("r.dat");
        }
        
        double ee = atoms.e_elec();
        double ep = atoms.e_pair();
        double ek = atoms.e_kin();
        double e_tot = ep + ek + ee;
        
        double e_gap = atoms.compute_Eg();
        
        fs << r * dt << '\t' << e_tot << '\t' << ee << '\t' << ep << '\t' << ek << '\t';
        fs << endl;
        
        if(r % n_bb == 0) {
            double bott_idx = atoms.compute_Bott_index();
            fs2 << bott_idx << endl;
        }
        
        if(r % n_save == 0) {
            
            mddata.basic_measurements(atoms);
            mddata.compute_g_r(atoms.pairs, atoms.numAtoms, atoms.volume, r_min, r_max, atoms.potential->dimension);
            
            average_MD_Data(accu_data, mddata, nx);
            nx++;
            
            accu_data.print_g_r("gr.dat");
            
            accu_data.print_data("a.dat", density);
            
            atoms.save_configuration("pos.dat");
        }
        squaremean += e_tot;
        meansquare += e_tot*e_tot;
    }
    fs.close();
    fs2.close();
}


void analyze_structures(int N_atoms, double density) {
    
    double filling_fraction = 0.5;      // filling fraction
    double mass = 1.0;                    // atom mass
    double gamma = 0.1;                // Langevin damping coefficient
    
    double l_box = pow(N_atoms / density, 1./3.);
    cout << "l_box = " << l_box << endl;
    
    double kT = 1.e-10;
    Simple_TB_system atoms(N_atoms, mass, filling_fraction, kT, gamma);
    atoms.set_BD_type(1);
    
    cout << "LD_tag = " << atoms.langevin_dym << endl;
    cout << "kT = " << kT << endl;
    cout << "density = " << density << endl;
    cout << "N_atom = " << N_atoms << endl;
        
    cout << "TI-TB parameters: " << endl;
    double M = -0.5;
    double lambda = 0; // 0.5;
    double t2 = 0.0;
    atoms.potential->set_parameters(M, lambda, t2);
    cout << atoms.potential->M << '\t' << atoms.potential->lambda << '\t' << atoms.potential->t2 << endl;
    
    
    double lt = l_box;
    double lz = l_box;
    vec3 bd_lo, bd_hi;
    
    //----- initialize the system ------
    atoms.init_random(lt, lz, bd_lo, bd_hi);
    //atoms.init_from_file(lt, lz, bd_lo, bd_hi);
    
    atoms.init_potential();

    atoms.find_all_pairs();
    
    int numOrbs = atoms.potential->numOrbitalsPerSite();
    ofstream fs("dlt.dat");
    for(auto p : atoms.pairs) {
        
        arma::cx_mat h(numOrbs, numOrbs);
        arma::cx_mat dhdx(numOrbs, numOrbs);
        arma::cx_mat dhdy(numOrbs, numOrbs);
        arma::cx_mat dhdz(numOrbs, numOrbs);
        
        atoms.potential->fill_TB_hoppings(p, h, dhdx, dhdy, dhdz);
        
        int ix = 2;
        if(p.index1 == ix || p.index2 == ix) {
            
            fs << p.delta.x << '\t' << p.delta.y << '\t' << p.delta.norm() << '\t' << abs(h(0,0)) << '\t';
            fs << h(0,0).real() << '\t' << h(1,1).real() << '\t';
            fs << h(0,1).real() << '\t' << h(0,1).imag() << '\t';
            fs << h(1,0).real() << '\t' << h(1,0).imag() << endl;
            
        }
    }
    fs.close();

    atoms.built_Hamiltonian();
    atoms.compute_density_matrix();
    atoms.print_eig_E("eig.dat");

    double e_g = atoms.compute_Eg();
    cout << "E_gap = " << e_g << endl;
    
    double bt_idx = atoms.compute_Bott_index();
    cout << "Bott idx = " << bt_idx << endl;
}

int main(int argc, const char * argv[]) {
    std::cout << "Hello, World!\n";

    //test_Bott_idx();
    
    double density          = argc > 1 ? atof(argv[1]) : 0.1;
    int N                   = argc > 2 ? atoi(argv[2]) : 300;
    double kB_T             = argc > 3 ? atof(argv[3]) : 0.05;
    double r_kT             = argc > 4 ? atof(argv[4]) : 0.02;
    
    cout << "kT_init = " << kB_T << "\t r_kT = " << r_kT << endl;
    
    int N_steps = 10000;
    double dt = 0.01;                   // time step

    cout.precision(12);
    //test_all();

    
    kB_T = 1.e-8;
    //QMD_simulation_noneq(N, density, kB_T, N_steps, dt);
    QMD_simulation(N, density, kB_T, N_steps, dt);

    //QMD_annealing_simulation(N, density, kB_T, r_kT, N_steps, dt);
    /*ofstream fs ("numvar.dat");
    for (int i = 1; i<=40; i++)
    {
        N=i;
        variance = 0; 
        for (int j = 0; j<20; j++) 
        {
            variance += QMD_simulation_noneq(N, density, kB_T, N_steps, dt);
        }
        variance = variance / 20.0; 
        fs << N << '\t' << variance << '\t';
        fs << endl;
    }
    fs.close();*/
    
    //QMD_simulation_open_BC(N, density, kB_T, N_steps, dt);

    //analyze_structures(200, 0.026250);
    
    return 0;
}

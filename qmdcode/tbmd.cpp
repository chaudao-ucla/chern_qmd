//
//  tbmd.cpp
//  TI-QMD
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#include <cassert>
#include "tbmd.hpp"

#include <armadillo>


//=========================================================
// Constructor
//=========================================================

Simple_TB_system::Simple_TB_system(int nAtoms, double atomMass, double filling, 
    double temperature, double gm) {
    numAtoms = nAtoms;
    position = vector<vec3>(nAtoms);
    velocity = vector<vec3>(nAtoms);
    force = vector<vec3>(nAtoms);
    
    mass = atomMass;
    kT = temperature;
    gamma = gm;
    
    ee = ep = ek = 0;
    
    kT_electron = kT;
    
    boundary_type = 0;
    langevin_dym = 1;
    
    
    //    Hamiltonian = arma::sp_cx_mat(dimH, dimH);
    //    rho = arma::sp_cx_mat(dimH, dimH);
    
    //potential = new(ToyModelSOribtal);
    potential = new(TI_2D_atom);
    
    dimH = nAtoms * potential->numOrbitalsPerSite();
    
    eigE = arma::vec(dimH);
    eig_mode = arma::cx_mat(dimH, dimH);

    
    filling_fraction = filling;
    num_filling = (int) (filling_fraction * dimH);

    U_mat = arma::cx_mat(num_filling, num_filling);
    W_mat = arma::cx_mat(num_filling, num_filling);
    
    Vx_mat = arma::cx_mat(dimH, dimH);
    Vy_mat = arma::cx_mat(dimH, dimH);
    
    rng = std::mt19937(seed());
    build_super_cells(0.5*potential->rcut());

    int N_c;
    int L_c;
    double c_size;
    Super_Cell *cell;
}

//=========================================================
// Initialization Functions
//=========================================================


void Simple_TB_system::init_random(double Lt, double Lz, vec3& bdsLo, vec3& bdsHi) {
    //RNG rg(0);
    std::uniform_real_distribution<double> rd;
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);
    
    for(int i=0; i<numAtoms; i++) {
        position[i] = vec3 {rd(rng) * Lt, rd(rng) * Lt, rd(rng) * Lz};
        
        if(potential->dimension == 2) position[i].z = 0;
    }
    
    system_box = bdsHi - bdsLo;
    cout << "system box = " << system_box << endl;
    if(potential->dimension == 3)
        volume = system_box.x * system_box.y * system_box.z;
    else if(potential->dimension == 2)
        volume = system_box.x * system_box.y;
    
    cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    double sgm_v = sqrt(kT / mass);
    std::normal_distribution<double> rn;    // default mean = 0, var = 1
    for(int i=0; i<numAtoms; i++) {
        velocity[i].x = sgm_v * rn(rng);
        velocity[i].y = sgm_v * rn(rng);
        velocity[i].z = sgm_v * rn(rng);
        
        if(potential->dimension == 2) velocity[i].z = 0;
    }
    
    pos_init = vector<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void Simple_TB_system::init_from_file(double Lt, double Lz, vec3 &bdsLo, vec3 &bdsHi) {
    //RNG rg(0);
    std::uniform_real_distribution<double> rd;
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);

    ifstream fs("pos_init.dat");
    
    for(int i=0; i<numAtoms; i++) {
        
        float x, y;
        
        fs >> x >> y;
        position[i].x = x;
        position[i].y = y;
        position[i].z = 0;
        
        //if(potential->dimension == 2) position[i].z = 0;
    }
    
    fs.close();
 
    /*
    ofstream fx("pp.dat");
    for(int i=0; i<numAtoms; i++) {
        fx << "(" << position[i].x << ", " << position[i].y << ")" << endl;
    }
    fx.close();
    */
    
    system_box = bdsHi - bdsLo;
    std::cout << "system box = " << system_box << endl;
    if(potential->dimension == 3)
        volume = system_box.x * system_box.y * system_box.z;
    else if(potential->dimension == 2)
        volume = system_box.x * system_box.y;
    
    std::cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = vec3 {0, 0, 0};
        force[i] = vec3 {0, 0, 0};
    }
    
    pos_init = vector<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void Simple_TB_system::init_from_file_vel(double Lt, double Lz, vec3 &bdsLo, vec3 &bdsHi) {
    //RNG rg(0);
    std::uniform_real_distribution<double> rd;
    
    bdsLo = {0, 0, 0};
    bdsHi = {Lt, Lt, Lz};
    position.resize(numAtoms);
    velocity.resize(numAtoms);

    ifstream fs("pos_init.dat");
    
    for(int i=0; i<numAtoms; i++) {
        
        float x, y, velx, vely;
        
        fs >> x >> y >> velx >> vely;
        position[i].x = x;
        position[i].y = y;
        position[i].z = 0;
        velocity[i].x = velx;
        velocity[i].y = vely;
        velocity[i].z = 0;   
        //if(potential->dimension == 2) position[i].z = 0;
    }
    fs.close();
 
    /*
    ofstream fx("pp.dat");
    for(int i=0; i<numAtoms; i++) {
        fx << "(" << position[i].x << ", " << position[i].y << ")" << endl;
    }
    fx.close();
    */
    
    system_box = bdsHi - bdsLo;
    std::cout << "system box = " << system_box << endl;
    if(potential->dimension == 3)
        volume = system_box.x * system_box.y * system_box.z;
    else if(potential->dimension == 2)
        volume = system_box.x * system_box.y;
    
    std::cout << "volumn = " << volume << endl;
    
    for(int i=0; i<numAtoms; i++) {
        force[i] = vec3 {0, 0, 0};
    }
    
    pos_init = vector<vec3>(numAtoms);
    for(int i=0; i<numAtoms; i++) pos_init[i] = position[i];
}

void Simple_TB_system::init_potential(void) {
    find_all_pairs();
    built_Hamiltonian();
    compute_density_matrix();
}

//=========================================================
// Pair Finding Algorithm
//=========================================================

void Simple_TB_system::find_all_pairs(void) {
    if(boundary_type == 0)
        pairs = allPairs(position, potential->rcut());
    else
        pairs = allPairsPeriodic(position, potential->rcut(), system_box);
}

vector<AtomPair> Simple_TB_system::get_pairs(vector<vec3> pos) {
    if(boundary_type == 0)
        return allPairs(pos, potential->rcut());
    else
        return allPairsPeriodic(pos, potential->rcut(), system_box);
}

void Simple_TB_system::find_all_pairs(double cutoff) {
    if(boundary_type == 0)
        pairs = allPairs(position, cutoff);
    else
        pairs = allPairsPeriodic(position, cutoff, system_box);
}

void Simple_TB_system::build_super_cells(double r_c) 
{
    
    L_c = (int) (system_box.x / r_c); // Assumes that the box is a square
    N_c = L_c * L_c;
    c_size = system_box.x / ((double) L_c);
    
    cout << "L_c = " << L_c << "\t N_c = " << N_c << "\t c_size = " << c_size << endl;
    
    cell = new Super_Cell[N_c];
    
    for(int x=0; x<L_c; x++)
    for(int y=0; y<L_c; y++)
    {
        int z = 0;
        int i = index_cell(x, y, z);
        
        cell[i].idx = i;
        cell[i].x = x;
        cell[i].y = y;
        cell[i].z = z;
        
        cell[i].r_min[0] = x * c_size;
        cell[i].r_max[0] = (x + 1.) * c_size;
        
        cell[i].r_min[1] = y * c_size;
        cell[i].r_max[1] = (y + 1.) * c_size;
        
        cell[i].r_min[2] = z * c_size;
        cell[i].r_max[2] = (z + 1.) * c_size;
        
        cell[i].atoms.clear();
    }
    
    for(int i=0; i<N_c; i++) {
        
        int x = cell[i].x;
        int y = cell[i].y;
        int z = cell[i].z;
        
        int j;
        
        int k = 0;
        
        int *tag = new int[N_c];
        for(int m=0; m<N_c; m++) tag[m] = 0;
        tag[i] = 1;
        
        // ====================================================
        
        j = index_cell(mod(x-1, L_c), y, z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), y, z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y-1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y+1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, y, mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, y, mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        // ====================================================
        
        j = index_cell(mod(x+1, L_c), mod(y+1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), mod(y-1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y+1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y-1, L_c), z);
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y+1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y+1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y-1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(x, mod(y-1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), y, mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), y, mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), y, mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), y, mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        // ====================================================
        
        j = index_cell(mod(x+1, L_c), mod(y+1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), mod(y-1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y+1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y-1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y-1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x-1, L_c), mod(y+1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), mod(y-1, L_c), mod(z+1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }
        
        j = index_cell(mod(x+1, L_c), mod(y+1, L_c), mod(z-1, L_c));
        if(tag[j] == 0) {
            cell[i].nn[k++] = &cell[j];
            tag[j] = 1;
        }

        delete [] tag;

        cell[i].N_neighbor = k;
    }
}

void Simple_TB_system::build_cell_atoms_list(void) 
{
    
    for(int r=0; r<N_c; r++) cell[r].atoms.clear();

    for(int i=0; i<Ns; i++) {
        double rx[3] = {position[i].x, position[i].y, position[i].z};
        int id[3];
        for(int m=0; m<3; m++) {
            id[m] = (int) (rx[m] / c_size);
        }
        
        int r = index_cell(id[0], id[1], id[2]);
        
        //        cout << "rx = " << rx[0] << '\t' << rx[1] << '\t' << rx[2] << endl;
        //        cout << cell[r].r_min[0] << '\t' << cell[r].r_min[1] << '\t' << cell[r].r_min[2] << endl;
        //        cout << cell[r].r_max[0] << '\t' << cell[r].r_max[1] << '\t' << cell[r].r_max[2] << endl;
                
        assert(rx[0] <= cell[r].r_max[0] && rx[0] >= cell[r].r_min[0] && rx[1] <= cell[r].r_max[1] && rx[1] >= cell[r].r_min[1] && rx[2] <= cell[r].r_max[2] && rx[2] >= cell[r].r_min[2]);
        
        cell[r].atoms.push_back(i);
        
        atom_index[i] = r;
    }
}

void Simple_TB_system::compute_neighbor_lists(void) {
    
    for(int i=0; i<Ns; i++) {
        neighbor_list[i].clear();
        
        int r = atom_index[i];
        
        int j;
        vec3 delta;
        double dr;
        
        for(int k=0; k<cell[r].atoms.size(); k++) {
            j = cell[r].atoms[k];
            if(j != i) {
                delta = wrapDelta(atom[j].position - atom[i].position, system_box);
                dr = delta.norm();
                
                if(dr < r_max) {
                    Neighbor tmp;
                    tmp.idx = j;
                    tmp.delta = delta;
                    tmp.dr = dr;
                    atom[i].neighbor_list.push_back(tmp);
                }
            }
        }
        
        for(int m=0; m<cell[r].N_neighbor; m++) {
            for(int k=0; k<cell[r].nn[m]->atoms.size(); k++) {
                j = cell[r].nn[m]->atoms[k];
                delta = wrapDelta(atom[j].position - atom[i].position, system_box);
                dr = delta.norm();
                
                if(dr < r_max) {
                    Neighbor tmp;
                    tmp.idx = j;
                    tmp.delta = delta;
                    tmp.dr = dr;
                    atom[i].neighbor_list.push_back(tmp);
                    
                }
            }
        }
        
    }
}

//=========================================================
// Build Hamiltonian and Computations
//=========================================================

void Simple_TB_system::built_Hamiltonian(void) {
    Hamiltonian = potential->build_Hamiltonian(numAtoms, pairs);
    //Hamiltonian = potential->buildHamiltonian(numAtoms, pairs, system_box);
}

double Simple_TB_system::compute_avg_density(double x) {
    double sum = 0;
    for(int i=0; i<eigE.size(); i++) {
        sum += fermi_density(eigE(i), kT, x);
    }
    return sum/((double) numAtoms);
}

double Simple_TB_system::compute_chemical_potential(void) {
    
    double x1 = eigE(0);
    double x2 = eigE(eigE.size()-1);
    
    int max_bisection = 100;
    double eps_bisection = 1.e-12;
    
    int iter = 0;
    while(iter < max_bisection || fabs(x2 - x1) > eps_bisection) {
        
        double xm = 0.5*(x1 + x2);
        double density = compute_avg_density(xm);
        
        if(density <= filling_fraction) x1 = xm;
        else x2 = xm;
        
        iter++;
    }
    
    return 0.5*(x1 + x2);
}
//=========================================================
// recomputes density matrix for each snapshot of time-independent simulation
//=========================================================
arma::cx_mat Simple_TB_system::compute_density_matrix(arma::cx_mat & H_elec) {
    
    arma::cx_mat rho_elec(H_elec.n_rows, H_elec.n_rows);
    fermi_energies.resize(H_elec.n_rows);
    arma::cx_mat eigvec;
    auto Hd = H_elec; //fkpm::sparseToDense(H_elec);
    arma::eig_sym(eigE, eigvec, Hd);
    eig_mode = eigvec;
    proj_states = eigvec;
    mu = compute_chemical_potential();
    
    for(int a=0; a<dimH; a++)
    for(int b=0; b<dimH; b++) {
        cx_double sum = 0;
        for(int i=0; i<eigE.size(); i++) {
            fermi_energies(i) = fermi_density(eigE(i), kT_electron, mu);
            sum += fermi_energies(i) * conj(eigvec(b, i)) * eigvec(a, i);
        }
        rho_elec(a, b) = sum;
    }
    
    return rho_elec;
}

//=========================================================
// first order Runge-Kutta 
//=========================================================

arma::cx_mat Simple_TB_system::evolve_density_matrix_RK1(arma::cx_mat & H_elec, 
        arma::cx_mat & proj_states, arma::cx_vec & fermi_energies, double dt)
{
    arma::cx_mat rho_elec(H_elec.n_rows, H_elec.n_rows);
    arma::cx_mat Hd = H_elec;  
    arma::eig_sym(eigE, eig_mode, Hd);
    arma::cx_mat update_proj(proj_states.n_rows,proj_states.n_cols);

    for (int i=0; i<proj_states.n_cols;i++)
    {
        arma::cx_mat temp(proj_states.n_rows,1);
        temp = proj_states.col(i);
        temp+= _I*-1.0*H_elec*dt;
        update_proj.col(i) = temp;
    }
    proj_states = update_proj;

    for(int a=0; a<dimH; a++)
    for(int b=0; b<dimH; b++) {
        cx_double sum = 0;
        for(int i=0; i<fermi_energies.size(); i++) {
            sum += fermi_energies(i) * conj(update_proj(b, i)) * update_proj(a, i);
        }
        rho_elec(a, b) = sum;
    }
    return rho_elec;
}

void Simple_TB_system::evolve_density_matrix_RK1(double dt)
{
    rho = evolve_density_matrix_RK1(Hamiltonian,proj_states,fermi_energies,dt);
}

void Simple_TB_system::compute_density_matrix(void) {
    rho = compute_density_matrix(Hamiltonian);
}

void Simple_TB_system::proj_to_rho(arma::cx_mat const & proj, arma::cx_mat & denmat)
{
    for(int a=0; a<dimH; a++)
    for(int b=0; b<dimH; b++) {
        cx_double sum = 0;
        for(int i=0; i<fermi_energies.size(); i++) {
            sum += fermi_energies(i) * conj(proj(b, i)) * proj(a, i);
        }
        denmat(a, b) = sum;
    }
}

void Simple_TB_system::compute_forces(void) {
    potential->force(rho, pairs, force, virial);
}

//void Simple_TB_system::compute_forces_elec(void) {
//    potential->force_elec(rho, pairs, force, virial);
//}

void Simple_TB_system::move_atoms(double dt) {
    for(int i=0; i<numAtoms; i++) {
        vec3 dlt = 0.5 * (force[i]/mass) * dt;
        dlt += velocity[i];
        velocity[i] = dlt;  // velecity at t + 0.5*dt
        dlt *= dt;
        position[i] += dlt;
    }
}
// Explicit first order Runge-Kutta to update positions
void Simple_TB_system::move_atoms_RK1(double dt)
{
    for (int i = 0; i<numAtoms; i++)
    {
        vec3 dlt = velocity[i]*dt;
        position[i]+=dlt;
    }
}

void Simple_TB_system::integrate_forces_RK1(double dt)
{
    for (int i=0; i<numAtoms; i++)
    {
        velocity[i]+=(force[i]/mass)*dt;
    }
}

void Simple_TB_system::integrate_forces(double dt) {
    for(int i=0; i<numAtoms; i++) {
        velocity[i] += 0.5 * (force[i]/mass) * dt;
    }
}

void Simple_TB_system::integrate_Langevin(double dt) {
    //RNG rng(seed());
    std::normal_distribution<double> rd;    // default mean = 0, var = 1
    
    double alpha2 = exp(-gamma * dt);
    double sigma2 = sqrt((1-pow(alpha2, 2))*kT/mass);
    for(int i=0; i<numAtoms; i++) {
        velocity[i] = alpha2 * velocity[i] + sigma2 * vec3(rd(rng), rd(rng), rd(rng));
        
        if(potential->dimension == 2) velocity[i].z = 0;
    }
}

void Simple_TB_system::step_NVT(double dt) {
    
    // velocity verlet integration with Langevin damping
    move_atoms(dt);
    find_all_pairs();
    
    built_Hamiltonian();
    compute_density_matrix();
    
    compute_forces();
    integrate_forces(dt);
    if(langevin_dym == 1) integrate_Langevin(dt);
}

void Simple_TB_system::evolve_NVT_RK1(double dt)
{
    move_atoms_RK1(dt);
    integrate_forces_RK1(dt);
    evolve_density_matrix_RK1(dt);

    //set-up for next round evolution step
    find_all_pairs();
    built_Hamiltonian();
    compute_forces();
}

//=========================================================
// Von-Neumann method of evolving density matrix 
//=========================================================
void Simple_TB_system::evolve_NVT_RK4_v2(double dt)
{
    arma::cx_mat H(Hamiltonian);
    arma::cx_mat P(proj_states);
    arma::cx_mat D(rho);
    arma::cx_mat D2, KD_sum, KD;
    arma::cx_mat KP_sum;
    arma::cx_mat P2, KP;

    vector<vec3> kforces(numAtoms);
    potential->force(rho,pairs,kforces,virial);

    vector<vec3> kpos(numAtoms);
    vector<vec3> kvel(numAtoms);

    vector<vec3> pos2(numAtoms);
    vector<vec3> vel2(numAtoms);

    vector<vec3> pos_sum(numAtoms);
    vector<vec3> vel_sum(numAtoms);

    // ------- RK4 step-1: ----------------

    KP = -_I * dt * H * P;
    P2 = P + 0.5 * KP;
    KP_sum = KP / 6.;

    KD = -_I * dt * (H * D - D * H);
    D2 = D + 0.5 * KD;
    KD_sum = KD / 6.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * velocity[i];
        pos_sum[i] = kpos[i] / 6.;
        pos2[i] = position[i] + 0.5 * kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] = kvel[i] / 6.;
        vel2[i] = velocity[i] + 0.5 * kvel[i];
    }

    // ------- RK4 step-2: ----------------

    vector<AtomPair> tempPairs = get_pairs(pos2);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(D2,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    P2 = P + 0.5 * KP;
    KP_sum += KP / 3.;

    KD = -_I * dt * (H * D2 - D2 * H);
    D2 = D + 0.5 * KD;
    KD_sum += KD / 3.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 3.;
        pos2[i] = position[i] + 0.5 * kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 3.;
        vel2[i] = velocity[i] + 0.5 * kvel[i];
    }

    // ------- RK4 step-3: ----------------

    tempPairs = get_pairs(pos2);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(D2,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    P2 = P + KP;
    KP_sum += KP / 3.;

    KD = -_I * dt * (H * D2 - D2 * H);
    D2 = D + KD;
    KD_sum += KD / 3.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 3.;
        pos2[i] = position[i] + kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 3.;
        vel2[i] = velocity[i] + kvel[i];
    }

    // ------- RK4 step-4: ----------------

    tempPairs = get_pairs(pos2);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(D2,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    KP_sum += KP / 6.;

    KD = -_I * dt * (H * D2 - D2 * H);
    KD_sum += KD / 6.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 6.;

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 6.;
    }

    // compute projected states, density matrix, Hamiltonian, position, velocity

    proj_states += KP_sum;

    for (int i=0; i<numAtoms; i++)
    {
        position[i] += pos_sum[i];
        velocity[i] += vel_sum[i];
    }

    rho = D + KD_sum;
    
    find_all_pairs();
    built_Hamiltonian();
}

//=========================================================
// Schrodinger equation method of evolving density matrix
//=========================================================

void Simple_TB_system::evolve_NVT_RK4(double dt)
{
    arma::cx_mat H(Hamiltonian);
    arma::cx_mat P(proj_states);
    arma::cx_mat KP_sum;
    arma::cx_mat P2, KP;

    vector<vec3> kforces(numAtoms);
    potential->force(rho,pairs,kforces,virial);

    vector<vec3> kpos(numAtoms);
    vector<vec3> kvel(numAtoms);

    vector<vec3> pos2(numAtoms);
    vector<vec3> vel2(numAtoms);

    vector<vec3> pos_sum(numAtoms);
    vector<vec3> vel_sum(numAtoms);

    // ------- RK4 step-1: ----------------

    KP = -_I * dt * H * P;
    P2 = P + 0.5 * KP;
    KP_sum = KP / 6.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * velocity[i];
        pos_sum[i] = kpos[i] / 6.;
        pos2[i] = position[i] + 0.5 * kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] = kvel[i] / 6.;
        vel2[i] = velocity[i] + 0.5 * kvel[i];
    }

    // ------- RK4 step-2: ----------------

    arma::cx_mat rho_temp(dimH,dimH);
    proj_to_rho(P2,rho_temp);

    vector<AtomPair> tempPairs = get_pairs(pos2);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(rho_temp,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    P2 = P + 0.5 * KP;
    KP_sum += KP / 3.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 3.;
        pos2[i] = position[i] + 0.5 * kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 3.;
        vel2[i] = velocity[i] + 0.5 * kvel[i];
    }

    // ------- RK4 step-3: ----------------

    tempPairs = get_pairs(pos2);
    proj_to_rho(P2,rho_temp);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(rho_temp,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    P2 = P + KP;
    KP_sum += KP / 3.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 3.;
        pos2[i] = position[i] + kpos[i];

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 3.;
        vel2[i] = velocity[i] + kvel[i];
    }

    // ------- RK4 step-4: ----------------

    tempPairs = get_pairs(pos2);
    proj_to_rho(P2,rho_temp);
    H = potential->build_Hamiltonian(numAtoms,tempPairs);
    potential->force(rho_temp,tempPairs,kforces,virial);

    KP = -_I * dt * H * P2;
    KP_sum += KP / 6.;

    for (int i=0; i<numAtoms; i++)
    {
        kpos[i] = dt * vel2[i];
        pos_sum[i] += kpos[i] / 6.;

        kvel[i] = (kforces[i]/mass)*dt;
        vel_sum[i] += kvel[i] / 6.;
    }

    // compute projected states, density matrix, Hamiltonian, position, velocity

    proj_states += KP_sum;

    for (int i=0; i<numAtoms; i++)
    {
        position[i] += pos_sum[i];
        velocity[i] += vel_sum[i];
    }

    proj_to_rho(proj_states,rho);
    find_all_pairs();
    built_Hamiltonian();
}

void Simple_TB_system::scale(double factor)
{
    vec3 bdsLo = {0, 0, 0};
    vec3 bdsHi = system_box;
    for(int i=0; i<numAtoms; i++) {
        position[i] = wrapPosition(position[i], bdsLo, bdsHi);
        position[i].x = position[i].x * factor;
        position[i].y = position[i].y * factor;
        position[i].z = position[i].z * factor;
    }
    system_box = bdsHi * factor;
    cout << "system box = " << system_box << endl;
    if(potential->dimension == 3)
        volume = system_box.x * system_box.y * system_box.z;
    else if(potential->dimension == 2)
        volume = system_box.x * system_box.y;
    cout << "volumn = " << volume << endl;

    find_all_pairs();
    built_Hamiltonian();
    compute_density_matrix();
}

double Simple_TB_system::e_elec(void) {
    double sum1 = 0;
    
    //for(int i=0; i<numFilling; i++) sum1 += eigE(i);
    
    //    cout << "mu = " << mu << endl;
    //
    //            for(int i=0; i<dimH; i++) {
    //                sum1 += eigE(i) * fermiDensity(eigE(i), kT, mu);
    //            }
    
    auto rr = Hamiltonian * rho;
    sum1 = trace(rr).real();
    
    ee = potential->numSpins() * sum1 / ((double) numAtoms);
    
    return ee;
}

double Simple_TB_system::e_kin(void) {
    double sum = 0;
    for(int i=0; i<numAtoms; i++) sum += velocity[i].norm2();
    ek = (0.5 * mass * sum) / ((double) numAtoms);
    return ek;
}

double Simple_TB_system::e_pair(void) {
    ep = potential->pair_energy(pairs) / ((double) numAtoms);
    return ep;
}

double Simple_TB_system::compute_Eg(void) {
    
    e_gap = eigE(num_filling) - eigE(num_filling - 1);
    return e_gap;
}

double Simple_TB_system::compute_Bott_index(void) {
    
    U_mat.zeros();
    W_mat.zeros();
    
    int numOrbs = potential->numOrbitalsPerSite();
    
    for(int a=0; a<num_filling; a++)
    for(int b=0; b<num_filling; b++) {
        cx_double sum1 = 0;
        cx_double sum2 = 0;
        
        for(int i=0; i<numAtoms; i++) {
            
            double theta = 2. * PI * (position[i].x / system_box.x);
            double phi   = 2. * PI * (position[i].y / system_box.y);
            
            for(int r=0; r<numOrbs; r++) {
                
                int m = i*numOrbs + r;
                
                sum1 += conj(proj_states(m, a)) * proj_states(m, b) * exp(_I * theta);
                sum2 += conj(proj_states(m, a)) * proj_states(m, b) * exp(_I * phi);
            }
        }
        
        U_mat(a, b) = sum1;
        W_mat(a, b) = sum2;
    }
    
    arma::cx_mat A_mat = W_mat * U_mat * W_mat.t() * U_mat.t();
    
    arma::cx_mat B_mat = arma::logmat(A_mat);
    
    
    cx_double tr = arma::trace(B_mat);
    
    //cout << "tr = " << tr << endl;
    
    
    return tr.imag() / (2. * PI);
}

void Simple_TB_system::compute_velocity_matrix(void) {
    
    Vx_mat.zeros();
    Vy_mat.zeros();
    
    double rc2 = pow(potential->rcut(), 2);
    
    int numOrbs = potential->numOrbitalsPerSite();

    for (auto const& p : pairs) {
        if (p.delta.norm2() < rc2) {
            int i = p.index1;
            int j = p.index2;
            assert(i < j);

            for(int o1 = 0; o1 < numOrbs; o1++)
            for(int o2 = 0; o2 < numOrbs; o2++) {

                int a = i*numOrbs+o1;
                int b = j*numOrbs+o2;
                
                Vx_mat(a, b) = -_I * p.delta.x * Hamiltonian(a, b);
                Vy_mat(a, b) = -_I * p.delta.y * Hamiltonian(a, b);
                
                Vx_mat(b, a) = +_I * p.delta.x * Hamiltonian(b, a);
                Vy_mat(b, a) = +_I * p.delta.y * Hamiltonian(b, a);
                    
            }
        }
    }
}

double _EPS_SIGMA = 1.e-5;

double Simple_TB_system::compute_sigma_xy(void) {
    
    compute_velocity_matrix();
    double sum = 0;
    
    for(int a=0; a<num_filling; a++) {
        for(int b=num_filling; b<dimH; b++) {
            
            
            cx_double v_x = arma::cdot(eig_mode.col(a), Vx_mat * eig_mode.col(b));
            cx_double v_y = arma::cdot(eig_mode.col(b), Vy_mat * eig_mode.col(a));
            
            sum += -2. * imag(v_x * v_y) / pow(eigE(a) - eigE(b) + _EPS_SIGMA, 2);
        }
    }
    
    return 2. * PI * sum / volume;
}

void Simple_TB_system::save_configuration(const string filename) {
    
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    vec3 bdsLo = {0, 0, 0};
    vec3 bdsHi = {system_box.x, system_box.y, system_box.z};

    for(int i=0; i<numAtoms; i++) {
        vec3 pos = position[i];
        vec3 vel = velocity[i];
        pos = wrapPosition(pos, bdsLo, bdsHi);
        fs << pos.x << '\t' << pos.y << '\t' << vel.x << '\t' << vel.y << '\t' << endl;
    }
    
    fs.close();
}

void Simple_TB_system::print_eig_E(const string filename) {

    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);
    
    for(int i=0; i<dimH; i++) {
        fs << eigE(i) << endl;
    }
    
    fs.close();
}


void Simple_TB_system::plot_eigen_mode(const string filename, int m) {
    
    auto U = eig_mode.col(m);
    
    std::ofstream fs;
    
    fs.open(filename.c_str(), ios::out);
    fs.precision(12);

    vec3 bdsLo = {0, 0, 0};
    vec3 bdsHi = {system_box.x, system_box.y, system_box.z};
    
    for(int i=0; i<numAtoms; i++) {
        vec3 pos = position[i];
        pos = wrapPosition(pos, bdsLo, bdsHi);
        fs << pos.x << '\t' << pos.y << '\t' << pow(abs(U(i)), 2) << endl;
    }
    
    fs.close();
    
}

void Simple_TB_system::test_1() {
    
    int m = (int) (rand1() * dimH);
    
    auto u = Hamiltonian * eig_mode.col(m);
    auto v = eigE(m) * eig_mode.col(m);
    
    auto dw = u - v;
    
    cout << "tst: " << arma::norm(dw) << '\t' << arma::norm(u) << '\t' << arma::norm(eig_mode.col(m)) << endl;
    
}





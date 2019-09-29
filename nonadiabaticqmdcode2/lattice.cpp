//
//  lattice.cpp
//  GvN_MD_Hubbard
//
//  Created by Gia-Wei Chern on 10/4/17.
//  Copyright © 2017 Gia-Wei Chern. All rights reserved.
//

#include <stdio.h>

#include <cassert>
#include "model.hpp"

using namespace std;

double TB_model::wrapDelta1d(double dy, double boxLength) {
    double dx = dy;
    while (dx > boxLength/2)
        dx -= boxLength;
    while (dx < -boxLength/2)
        dx += boxLength;
    return dx;
}

vec3 TB_model::wrapDelta(vec3 delta, vec3 boxLengths) {
    vec3 dlt;
    dlt.x = wrapDelta1d(delta.x, boxLengths.x);
    dlt.y = wrapDelta1d(delta.y, boxLengths.y);
    dlt.z = wrapDelta1d(delta.z, boxLengths.z);
    return dlt;
}

double TB_model::wrapPosition1d(double y, double bdsLo, double bdsHi) {
    double length = bdsHi - bdsLo;
    double x = y;
    while (x > bdsHi) {
        x -= length;
    }
    while (x < bdsLo) {
        x += length;
    }
    return x;
}

vec3 TB_model::wrapPosition(vec3 p, vec3 bdsLo, vec3 bdsHi) {
    vec3 r;
    r.x = wrapPosition1d(p.x, bdsLo.x, bdsHi.x);
    r.y = wrapPosition1d(p.y, bdsLo.y, bdsHi.y);
    r.z = wrapPosition1d(p.z, bdsLo.z, bdsHi.z);
    return r;
}

void TB_model::move_atoms_pbc(void) {
    vec3 bdsLo, bdsHi;
    bdsLo = {0, 0, 0};
    bdsHi = {system_box.x, system_box.y, system_box.z};
    
    for(int i=0; i<Ns; i++) {
        atom[i].position =  wrapPosition(atom[i].position, bdsLo, bdsHi);
    }
}

void TB_model::build_super_cells(double r_c) {
    
    L_c = (int) (L / r_c);
    N_c = L_c * L_c * L_c;
    c_size = L / ((double) L_c);
    
    cout << "L_c = " << L_c << "\t N_c = " << N_c << "\t c_size = " << c_size << endl;
    
    cell = new Super_Cell[N_c];
    
    for(int x=0; x<L_c; x++)
    for(int y=0; y<L_c; y++)
    for(int z=0; z<L_c; z++) {
                
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
    
//    ofstream fs("tst.dat");
//    for(int i=0; i<N_c; i++) {
//        fs << i << "(" << cell[i].x << ", " << cell[i].y << ", " << cell[i].z << ")" << endl;
//        fs << "----------------------------" << endl;
//        for(int m=0; m<cell[i].N_neighbor; m++) {
//            int j = cell[i].nn[m]->idx;
//            fs << j << "(" << cell[j].x << ", " << cell[j].y << ", " << cell[j].z << ")" << endl;
//        }
//        fs << "=====================================================" << endl << endl;
//    }
//    fs.close();
//    exit(1);
}

void TB_model::build_cell_atoms_list(void) {
    
    for(int r=0; r<N_c; r++) cell[r].atoms.clear();

    
    for(int i=0; i<Ns; i++) {
        double rx[3] = {atom[i].position.x, atom[i].position.y, atom[i].position.z};
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
        
        atom[i].cell_idx = r;
    }
}

void TB_model::compute_neighbor_lists(void) {
    
    for(int i=0; i<Ns; i++) {
        atom[i].neighbor_list.clear();
        
        int r = atom[i].cell_idx;
        
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

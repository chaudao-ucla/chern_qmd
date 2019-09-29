//
//  units.hpp
//  GQMD-Hubbard-liquid
//
//  Created by Gia-Wei Chern on 1/12/18.
//  Copyright © 2018 Gia-Wei Chern. All rights reserved.
//

#ifndef units_h
#define units_h

#include <armadillo>

// length
constexpr double angstrom = 1.0;
constexpr double bohr = 0.5292*angstrom; // 0.5291772109217

// energy
constexpr double eV = 1.0;
constexpr double rydberg = 13.605804*eV; // 13.6056925330

// mass
constexpr double amu = 1.0;
constexpr double massSi = 28.0851*amu;

// time
constexpr double femtosecond = /* sqrt(amu/eV) * angstrom */ 1. / 10.1805054836529;

// temperature
constexpr double kelvin = 1.0;
constexpr double kB = 8.6173323849609e-5 * eV / kelvin;

#endif /* units_h */

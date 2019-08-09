/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#pragma once

#include "mrchem.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/qmfunction_fwd.h"
#include <map>

namespace mrchem {

struct deposit {
    Orbital *orb;
    FunctionData funcinfo; // Function info
    OrbitalData orbinfo;   // Orbital info
    int id;                // to identify what is deposited
    int source;            // mpi rank from the source of the data
};

struct queue_struct {
    int id;
    std::vector<int> clients;
};

class Bank {
public:
    Bank();
    ~Bank();
    int open();
    int close();
    int clear_all(int i, MPI_Comm comm);
    int clear(int ix);
    int put_orb(int id, Orbital &orb);
    int get_orb(int id, Orbital &orb);
    int put_func(int id, QMFunction &func);
    int get_func(int id, QMFunction &func);

private:
    int const CLOSE_BANK = 1;
    int const CLEAR_BANK = 2;
    int const GET_ORBITAL = 3;
    int const SAVE_ORBITAL = 4;
    int const GET_FUNCTION = 5;
    int const SAVE_FUNCTION = 6;
    int clear_bank();
    std::map<int, int> id2ix;
    std::vector<deposit> deposits;
    std::map<int, int> id2qu;
    std::vector<queue_struct> queue;
};

} // namespace mrchem

/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <memory>
#include <MRCPP/MWOperators>
#include <MRCPP/trees/FunctionNode.h>
#include <XCFun/xcfun.h>
#include <xc_funcs.h>
#include <xc.h>

namespace mrdft {

using XC_p = std::unique_ptr<xcfun_t, decltype(&xcfun_delete)>;

class XClib {
public:
    xcfun_t *xcfun;                               ///< @brief XCFun library handle
    bool libxc;                                   ///< @brief Flag indicating if Libxc is active (True if "DFT {xc_library = libxc}" in input file)
    std::vector<xc_func_type*> libxc_objects;     ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;              ///< @brief Vector scaling coefficients for each functional in libxc_objects

};

} // mrdft
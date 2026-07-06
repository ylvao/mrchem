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

class XClib {
public:
    explicit XClib(bool spin_enabled) {
        this->spin = spin_enabled;
    }

    static bool libxc;                            ///< @brief Flag indicating if Libxc is active (True if "DFT {xc_library = libxc}" in input file)

    xcfun_t *xcfun;                               ///< @brief XCFun library handle
    std::vector<xc_func_type*> libxc_objects;     ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;              ///< @brief Vector scaling coefficients for each functional in libxc_objects
    bool spin;

    bool use_gamma_derivatives{false}; // from old code: remove???

    virtual double setFunctional(const std::string &name, double c, double cutoff) = 0;
    virtual double getAmountExx() const = 0;
    virtual int getnOut() = 0;
    virtual void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) = 0;
    virtual void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, double cutoff) const = 0;
    virtual void printFunctionalReference(int out_txt_width, std::vector<std::string> xcfun_func_names) const = 0;

    // Common printing helpers shared by Libxc and XCFun implementations
    static void printReferenceHeader(int out_txt_width);
    static void printWrap(const std::string &str, std::size_t txt_width, int indent = 0);

};

} // mrdft
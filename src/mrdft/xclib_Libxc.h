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

#include <xc_funcs.h>
#include <xc.h>
#include "xclib.h"

namespace mrdft {

/**
 * @class Libxc
 * @brief Class for Exchange-Correlation functionals using LibXC
 * @details Provides the interface for evaluating XC potentials
 * on the Multi-Resolution Analysis (MRA) grid using Libxc
 */
class Libxc final : public XClib {

public:
    explicit Libxc(bool spin_enabled)
        : XClib(spin_enabled) {}

    ~Libxc() override;

    std::vector<xc_func_type*> libxc_objects;     ///< @brief Vector of initialized Libxc functionals
    std::vector<double> libxc_coefs;              ///< @brief Vector scaling coefficients for each functional in libxc_objects
    double customExx{0.0};                        ///< @brief Used in mapfunctionalName to set exx for custom functionals
    void setFunctional(const std::string &name, double c, double cutoff) override;
    void setCustomExx(double exx) { customExx = exx; }
    double getAmountExx() const override;
    double getCustomExx() const { return customExx; }
    int getnOut() override;

    void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) override;
    void printFunctionalReference(int out_txt_width) const override;
    void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, double cutoff) const override;
    void mapFunctionalName(std::string name, std::vector<int> &ids, std::vector<double> &coefs, double &customExx);

};

} // mrdft
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
#include "xclib.h"

namespace mrdft {

/**
 * @class XCFun
 * @brief Class for Exchange-Correlation functionals using XCFun
 * @details Provides the interface for evaluating XC potentials
 * on the Multi-Resolution Analysis (MRA) grid using XCFun
 */
class XCFun final : public XClib {

public:
    double setFunctional(const std::string &name, double c, double cutoff, bool spin) override;
    double getAmountExx() const override;
    void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, bool spin, int order, bool gamma) override;
    void printFunctionalReference(int out_txt_width, std::vector<std::string> xcfun_func_names) const override;
    void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, bool spin, double cutoff) const override; 

};

} // mrdft
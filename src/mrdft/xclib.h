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

#include <MRCPP/MWOperators>

namespace mrdft {

class XClib {
public:
    explicit XClib(bool spin_enabled) { this->spin = spin_enabled; }

    virtual ~XClib() = default;

    std::vector<std::pair<std::string, double>> functional_specs; ///< @brief Configured functionals and scaling coefficients

    virtual void setFunctional(const std::string &name, double c, double cutoff) = 0;
    virtual double getAmountExx() const = 0;
    virtual int getnOut() = 0;
    virtual void initFunctionalLibrary(bool &lda, bool &gga, bool &mgga, int order, bool gamma) = 0;
    virtual void callLibEval(const Eigen::MatrixXd &inp, Eigen::MatrixXd &out, int nPts, int nInp, int nOut, double cutoff) const = 0;
    virtual void printFunctionalReference(int out_txt_width) const = 0;

    void addFunctionalSpec(const std::string &name, double c) { functional_specs.emplace_back(name, c); }

    static bool libxc; ///< @brief Flag indicating if Libxc is active (True if "DFT {xc_library = libxc}" in input file)
    bool spin;
    const std::vector<std::pair<std::string, double>> &getFunctionalSpecs() const { return functional_specs; }

    static void printReferenceHeader(int out_txt_width);
    static void printWrap(const std::string &str, std::size_t txt_width, int indent = 0);
};

} // namespace mrdft
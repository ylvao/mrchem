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
#include <memory> // added for use of libxc
#include <string> // added for use of libxc
#include <map> // added for use of libxc
#include <vector> // added for use of libxc
#include <xc.h> // added for use of libxc
// #include <XCFun/xcfun.h>

#include "MRDFT.h"

namespace mrdft {

// Structure to hold LibXC functional data
struct LibXCData {
    int func_id;             // Functional ID in LibXC
    double weight;           // Weight for this functional
    xc_func_type func;       // LibXC function type
    bool initialized;        // If the functional is initialized
};

class Factory final {
public:
    Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA);
    // ~Factory() = default;
    ~Factory()

    void setSpin(bool s) { spin = s; }
    void setOrder(int k) { order = k; }
    void setUseGamma(bool g) { gamma = g; }
    void setLogGradient(bool lg) { log_grad = lg; }
    void setDensityCutoff(double c) { cutoff = c; }
    void setDerivative(const std::string &n) { diff_s = n; }
    // void setFunctional(const std::string &n, double c = 1.0) { xcfun_set(xcfun_p.get(), n.c_str(), c); }

    // added for use of libxc
    void setFunctional(const std::string &name, double weight = 1.0);
    bool isGGA() const;
    bool isHybrid() const;
    bool getHybridCoeff() const;

    std::unique_ptr<MRDFT> build();

private:
    int order{1};
    bool spin{false};
    bool gamma{false};
    bool log_grad{false};
    double cutoff{-1.0};
    std::string diff_s{"abgv_00"};
    const mrcpp::MultiResolutionAnalysis<3> mra;

    // added for use of libxc
    std::vector<LibXCData> functionals;
    int mapFunctionalName(const std::string &name) const;
    void cleanupFunctionals();

    // XC_p xcfun_p;
    std::unique_ptr<mrcpp::DerivativeOperator<3>> diff_p;
};

} // namespace mrdft

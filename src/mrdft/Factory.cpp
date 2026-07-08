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

#include "Factory.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>

#include "GGA.h"
#include "Grid.h"
#include "LDA.h"
#include "MRDFT.h"
#include "SpinGGA.h"
#include "SpinLDA.h"
#include "xclib_Libxc.h"
#include "xclib_XCFun.h"

namespace mrdft {

bool XClib::libxc = false;

Factory::Factory(const mrcpp::MultiResolutionAnalysis<3> &MRA)
        : mra(MRA) {
    if (XClib::libxc) {
        xclib = std::make_unique<Libxc>(spin);
    } else {
        xclib = std::make_unique<XCFun>(spin);
    }
}
void Factory::rebuildXClib() {
    if (!xclib) {
        return;
    }

    const auto functional_specs = xclib->getFunctionalSpecs();
    if (XClib::libxc) {
        xclib = std::make_unique<Libxc>(spin);
    } else {
        xclib = std::make_unique<XCFun>(spin);
    }

    for (const auto &[name, coef] : functional_specs) {
        xclib->setFunctional(name, coef, cutoff);
    }
}

void Factory::setLibxc(bool libxc_) {
    if (XClib::libxc == libxc_) {
        return;
    }

    XClib::libxc = libxc_;
    rebuildXClib();
}

void Factory::setSpin(bool s) {
    if (spin == s) {
        return;
    }

    spin = s;
    rebuildXClib();
}


void Factory::setFunctional(const std::string &name, double c) {
    xclib->setFunctional(name, c, cutoff);
}

std::unique_ptr<MRDFT> Factory::build() {
    // Init DFT grid
    auto grid_p = std::make_unique<Grid>(mra);

    // Init Libxc or XCFun and set functional family bools
    bool gga = false;
    bool lda = false;
    bool mgga = false;
    xclib->initFunctionalLibrary(lda, gga, mgga, order, gamma);

    // Init MW derivative
    if (gga) {
        if (diff_s == "bspline") diff_p = std::make_unique<mrcpp::BSOperator<3>>(mra, 1);
        if (diff_s == "abgv_00") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
        if (diff_s == "abgv_55") diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.5, 0.5);
    }

    // Init XC functional
    std::unique_ptr<Functional> func_p{nullptr};
    if (spin) {
        if (gga)
            func_p = std::make_unique<SpinGGA>(order, xclib, diff_p);
        else if (lda)
            func_p = std::make_unique<SpinLDA>(order, xclib);
        else
            MSG_ABORT("Case not handled");
    } else {
        if (gga)
            func_p = std::make_unique<GGA>(order, xclib, diff_p);
        else if (lda)
            func_p = std::make_unique<LDA>(order, xclib);
        else
            MSG_ABORT("Case not handled");
    }
    if (func_p == nullptr) MSG_ABORT("Invalid functional type");
    diff_p = std::make_unique<mrcpp::ABGVOperator<3>>(mra, 0.0, 0.0);
    func_p->setDerivOp(diff_p);
    func_p->setLogGradient(log_grad);
    func_p->setDensityCutoff(cutoff);

    auto mrdft_p = std::make_unique<MRDFT>(grid_p, func_p);
    return mrdft_p;
}

} // namespace mrdft

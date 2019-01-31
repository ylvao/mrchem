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

namespace mrchem {

// clang-format off
class HyperFineCoupling final {
public:
    explicit HyperFineCoupling(const Nucleus &n) : nuc(n) {}

    const Nucleus &getNucleus() const { return this->nuc; }

    DoubleMatrix getTensor() const { return getSpinTerm() + getFermiContactTerm(); }
    DoubleMatrix &getSpinTerm() { return this->spin_term; }
    DoubleMatrix &getFermiContactTerm() { return this->fc_term; }
    const DoubleMatrix &getSpinTerm() const{ return this->spin_term; }
    const DoubleMatrix &getFermiContactTerm() const { return this->fc_term; }

    friend std::ostream& operator<<(std::ostream &o, const HyperFineCoupling &hfc) {
        auto fc_term = hfc.getFermiContactTerm()(0,0);
        auto spin_term = 1.0/hfc.getSpinTerm()(0,0);

        auto beta_e = PHYSCONST::beta_e;  // Bohr magneton
        auto beta_N = PHYSCONST::beta_N;  // Nuclear magneton
        auto g_e    = PHYSCONST::g_e;     // Free-electron g-value
        auto g_N    = 0.0; //hfc.getNucleus().getElement().getGValue();

        auto hfcc_g = 0.0;
        auto hfcc_hz = 0.0;

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                    HyperFine Coupling Constant             " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(5);
        o << std::setw(3)  << hfc.getNucleus().getElement().getSymbol();
        o << std::setw(26) << hfc.getNucleus().getCoord()[0];
        o << std::setw(15) << hfc.getNucleus().getCoord()[1];
        o << std::setw(15) << hfc.getNucleus().getCoord()[2];
        o << std::endl;
        mrcpp::Printer::setPrecision(10);
        o << "                                                            " << std::endl;
        o << "-------------------- Isotropic averages --------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " A                    (gauss)" << std::setw(30) << hfcc_g     << std::endl;
        o << "                      (MHz)  " << std::setw(30) << hfcc_hz    << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    const Nucleus nuc;
    DoubleMatrix fc_term{DoubleMatrix::Zero(1,1)};
    DoubleMatrix spin_term{DoubleMatrix::Zero(1,1)};
};
// clang-format on

} // namespace mrchem

/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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
class NMRShielding final {
public:
    explicit NMRShielding(const Nucleus &n) : nuc(n) {}

    const Nucleus& getNucleus() const { return this->nuc; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const NMRShielding &nmr) {
        auto isoDSppm = nmr.getDiamagnetic().trace() / 3.0;
        auto isoPSppm = nmr.getParamagnetic().trace() / 3.0;
        auto isoTSppm = isoDSppm + isoPSppm;

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                    NMR shielding tensor                    " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(5);
        o << std::setw(3)  << nmr.getNucleus().getElement().getSymbol();
        o << std::setw(26) << nmr.getNucleus().getCoord()[0];
        o << std::setw(15) << nmr.getNucleus().getCoord()[1];
        o << std::setw(15) << nmr.getNucleus().getCoord()[2];
        o << std::endl;
        mrcpp::Printer::setPrecision(10);
        o << "                                                            " << std::endl;
        o << "-------------------- Isotropic averages --------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Total            (ppm)      " << std::setw(30) << isoTSppm   << std::endl;
        o << " Diamagnetic      (ppm)      " << std::setw(30) << isoDSppm   << std::endl;
        o << " Paramagnetic     (ppm)      " << std::setw(30) << isoPSppm   << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Total ---------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                nmr.getTensor().format(clean_format)            << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Diamagnetic ------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                nmr.getDiamagnetic().format(clean_format)       << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Paramagnetic -----------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                nmr.getParamagnetic().format(clean_format)      << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    const Nucleus nuc;
    DoubleMatrix dia_tensor{DoubleMatrix::Zero(3,3)};
    DoubleMatrix para_tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} //namespace mrchem

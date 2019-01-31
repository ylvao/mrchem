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
class SpinSpinCoupling final {
public:
    SpinSpinCoupling(const Nucleus &n_k, const Nucleus &n_l) : nuc_K(n_k), nuc_L(n_l) {}

    const Nucleus &getNucleusK() const { return this->nuc_K; }
    const Nucleus &getNucleusL() const { return this->nuc_L; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const SpinSpinCoupling &sscc) {
        auto isoDShz = sscc.getDiamagnetic().trace() / 3.0;
        auto isoPShz = sscc.getParamagnetic().trace() / 3.0;
        auto isoTShz = isoDShz + isoPShz;

        int oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                    Spin-Spin Coupling Constant             " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(5);
        o << std::setw(3)  << sscc.getNucleusK().getElement().getSymbol();
        o << std::setw(26) << sscc.getNucleusK().getCoord()[0];
        o << std::setw(15) << sscc.getNucleusK().getCoord()[1];
        o << std::setw(15) << sscc.getNucleusK().getCoord()[2];
        o << std::endl;
        o << std::setw(3)  << sscc.getNucleusL().getElement().getSymbol();
        o << std::setw(26) << sscc.getNucleusL().getCoord()[0];
        o << std::setw(15) << sscc.getNucleusL().getCoord()[1];
        o << std::setw(15) << sscc.getNucleusL().getCoord()[2];
        o << std::endl;
        mrcpp::Printer::setPrecision(10);
        o << "                                                            " << std::endl;
        o << "-------------------- Isotropic averages --------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Total            (Hz)       " << std::setw(30) << isoTShz    << std::endl;
        o << " Diamagnetic      (Hz)       " << std::setw(30) << isoDShz    << std::endl;
        o << " Paramagnetic     (Hz)       " << std::setw(30) << isoPShz    << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Total ---------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                sscc.getTensor().format(clean_format)           << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Diamagnetic ------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                sscc.getDiamagnetic().format(clean_format)      << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Paramagnetic -----------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                sscc.getParamagnetic().format(clean_format)     << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    const Nucleus nuc_K;
    const Nucleus nuc_L;
    DoubleMatrix dia_tensor{DoubleMatrix::Zero(3,3)};
    DoubleMatrix para_tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} // namespace mrchem

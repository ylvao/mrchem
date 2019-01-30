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
class Magnetizability final {
public:
    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Magnetizability &mag) {
        auto w_au = 0.0;  // Only static magnetizability
        auto isoDMau = mag.getDiamagnetic().trace() / 3.0;
        auto isoPMau = mag.getParamagnetic().trace() / 3.0;
        auto isoTMau = isoDMau + isoPMau;

        auto isoDMsi = isoDMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})
        auto isoPMsi = isoPMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})
        auto isoTMsi = isoTMau * PHYSCONST::JT_m2; // SI units (J/T^2 10^{-30})

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                   Magnetizability tensor                   " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Frequency:       (au)       " << std::setw(30) << w_au       << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------- Isotropic averages --------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Total            (au)       " << std::setw(30) << isoTMau    << std::endl;
        o << " Diamagnetic      (au)       " << std::setw(30) << isoDMau    << std::endl;
        o << " Paramagnetic     (au)       " << std::setw(30) << isoPMau    << std::endl;
        o << "                                                            " << std::endl;
        o << " Total            (SI)       " << std::setw(30) << isoTMsi    << std::endl;
        o << " Diamagnetic      (SI)       " << std::setw(30) << isoDMsi    << std::endl;
        o << " Paramagnetic     (SI)       " << std::setw(30) << isoPMsi    << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Total ---------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                mag.getTensor().format(clean_format)            << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Diamagnetic ------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                mag.getDiamagnetic().format(clean_format)       << std::endl;
        o << "                                                            " << std::endl;
        o << "----------------------- Paramagnetic -----------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<                mag.getParamagnetic().format(clean_format)      << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    DoubleMatrix dia_tensor{DoubleMatrix::Zero(3,3)};
    DoubleMatrix para_tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} //namespace mrchem

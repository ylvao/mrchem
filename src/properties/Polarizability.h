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
class Polarizability final {
public:
    explicit Polarizability(double w = 0.0) : frequency(w) {}

    double getFrequency() const { return this->frequency; }
    DoubleMatrix &getTensor() { return this->tensor; }
    const DoubleMatrix &getTensor() const { return this->tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Polarizability &pol) {
        auto w_au = pol.getFrequency();
        auto iso_au = pol.getTensor().trace() / 3.0;
        auto iso_si = iso_au * 0.0; // Luca: FIX THIS

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                   Polarizability tensor                    " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Frequency:        (au)      " << std::setw(30) << w_au       << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Isotropic average (au)      " << std::setw(30) << iso_au     << std::endl;
        o << " Isotropic average (SI)      TO BE FIXED                    " << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Tensor --------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<             pol.getTensor().format(clean_format)               << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    double frequency;
    DoubleMatrix tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} //namespace mrchem

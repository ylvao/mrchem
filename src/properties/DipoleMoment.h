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
class DipoleMoment final {
public:
    DoubleVector getTensor() const { return getNuclear() + getElectronic(); }
    DoubleVector &getNuclear() { return this->nuc_tensor; }
    DoubleVector &getElectronic() { return this->el_tensor; }
    const DoubleVector &getNuclear() const { return this->nuc_tensor; }
    const DoubleVector &getElectronic() const { return this->el_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dip) {
        auto length_au = dip.getTensor().norm();
        auto length_db = length_au * PHYSCONST::Debye;

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                         Dipole moment                      " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Length of vector:   (au)    " << std::setw(30) << length_au  << std::endl;
        o << "                     (Debye) " << std::setw(30) << length_db  << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Total ---------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<    dip.getTensor().transpose().format(clean_format)            << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------- Nuclear --------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<    dip.getNuclear().transpose().format(clean_format)           << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------ Electronic ------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<    dip.getElectronic().transpose().format(clean_format)        << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    DoubleVector nuc_tensor{DoubleVector::Zero(3)};
    DoubleVector el_tensor{DoubleVector::Zero(3)};
};
// clang-format on

} // namespace mrchem

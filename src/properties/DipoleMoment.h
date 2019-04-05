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
#include "utils/math_utils.h"

namespace mrchem {

// clang-format off
class DipoleMoment final {
public:
    mrcpp::Coord<3> &getOrigin() { return this->origin; }
    const mrcpp::Coord<3> &getOrigin() const { return this->origin; }

    DoubleVector getTensor() const { return getNuclear() + getElectronic(); }
    DoubleVector &getNuclear() { return this->nuc_tensor; }
    DoubleVector &getElectronic() { return this->el_tensor; }
    const DoubleVector &getNuclear() const { return this->nuc_tensor; }
    const DoubleVector &getElectronic() const { return this->el_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const DipoleMoment &dip) {
        auto prec = mrcpp::Printer::getPrecision();
        auto length_au = dip.getTensor().norm();
        auto length_db = length_au * PHYSCONST::Debye;

        std::string origin_str = math_utils::coord_to_string(prec, 14, dip.getOrigin());

        std::stringstream o_au, o_db;
        o_au << std::setw(27) << std::setprecision(prec) << std::fixed << length_au;
        o_db << std::setw(27) << std::setprecision(prec) << std::fixed << length_db;

        std::string el_str = math_utils::vector_to_string(prec, 14, dip.getElectronic());
        std::string nuc_str = math_utils::vector_to_string(prec, 14, dip.getNuclear());
        std::string tot_str = math_utils::vector_to_string(prec, 14, dip.getTensor());

        o << "============================================================" << std::endl;
        o << "                         Dipole moment                      " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "            r_O :" << origin_str                              << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "     Electronic :" << el_str                                  << std::endl;
        o << "        Nuclear :" << nuc_str                                 << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Total vector :" << tot_str                                 << std::endl;
        o << "      Magnitude :          (au) " << o_au.str()               << std::endl;
        o << "                :       (Debye) " << o_db.str()               << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;

        return o;
    }

private:
    mrcpp::Coord<3> origin{};
    DoubleVector nuc_tensor{DoubleVector::Zero(3)};
    DoubleVector el_tensor{DoubleVector::Zero(3)};
};
// clang-format on

} // namespace mrchem

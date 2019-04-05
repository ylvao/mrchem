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
class Polarizability final {
public:
    explicit Polarizability(double w = 0.0) : frequency(w) {}

    mrcpp::Coord<3> &getOrigin() { return this->origin; }
    const mrcpp::Coord<3> &getOrigin() const { return this->origin; }

    double getFrequency() const { return this->frequency; }
    DoubleMatrix &getTensor() { return this->tensor; }
    const DoubleMatrix &getTensor() const { return this->tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Polarizability &pol) {
        auto prec = mrcpp::Printer::getPrecision();
        auto iso_au = pol.getTensor().trace() / 3.0;
        auto iso_si = iso_au * 0.0; // Luca: FIX THIS

        std::string origin_str = math_utils::coord_to_string(prec, 14, pol.getOrigin());

        std::stringstream o_omega;
        o_omega << std::setw(27) << std::setprecision(prec) << std::fixed << pol.getFrequency();

        std::stringstream o_iso_au, o_iso_si;
        o_iso_au << std::setw(27) << std::setprecision(prec) << std::fixed << iso_au;
        o_iso_si << std::setw(27) << std::setprecision(prec) << std::fixed << iso_si;

        std::string pol_str_0 = math_utils::vector_to_string(prec, 14, pol.getTensor().row(0));
        std::string pol_str_1 = math_utils::vector_to_string(prec, 14, pol.getTensor().row(1));
        std::string pol_str_2 = math_utils::vector_to_string(prec, 14, pol.getTensor().row(2));

        o << "============================================================" << std::endl;
        o << "                      Polarizability                        " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "      Frequency :          (au) " << o_omega.str()            << std::endl;
        o << "            r_O :" << origin_str                              << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Total tensor :" << pol_str_0                               << std::endl;
        o << "                :" << pol_str_1                               << std::endl;
        o << "                :" << pol_str_2                               << std::endl;
        o << "   Iso. average :          (au) " << o_iso_au.str()           << std::endl;
        o << "                :          (SI)                 TO BE FIXED " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;

        return o;
    }

private:
    double frequency;
    mrcpp::Coord<3> origin{};
    DoubleMatrix tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} // namespace mrchem

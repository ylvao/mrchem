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
#include "utils/print_utils.h"

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

    void print() const {
        auto iso_au = getTensor().trace() / 3.0;
        auto iso_si = iso_au * 0.0; // Luca: FIX THIS

        auto prec = mrcpp::Printer::getPrecision();
        mrcpp::print::header(0, "Polarizability");
        print_utils::scalar(0, "Frequency", "(au)", getFrequency(), prec, false);
        print_utils::coord(0, "Origin", getOrigin(), prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor(), prec, false);
        print_utils::scalar(0, "Iso. average", "(au)", iso_au, prec, false);
        print_utils::scalar(0, "            ", "(SI)", iso_si, prec, false);
        mrcpp::print::separator(0, '=', 2);
    }

private:
    double frequency;
    mrcpp::Coord<3> origin{};
    DoubleMatrix tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} // namespace mrchem

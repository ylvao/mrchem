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
class DipoleMoment final {
public:
    mrcpp::Coord<3> &getOrigin() { return this->origin; }
    const mrcpp::Coord<3> &getOrigin() const { return this->origin; }

    DoubleVector getTensor() const { return getNuclear() + getElectronic(); }
    DoubleVector &getNuclear() { return this->nuc_tensor; }
    DoubleVector &getElectronic() { return this->el_tensor; }
    const DoubleVector &getNuclear() const { return this->nuc_tensor; }
    const DoubleVector &getElectronic() const { return this->el_tensor; }

    void print() const {
        auto length_au = getTensor().norm();
        auto length_db = length_au * PHYSCONST::Debye;

        auto prec = mrcpp::Printer::getPrecision();
        mrcpp::print::header(0, "Dipole Moment");
        print_utils::coord(0, "Origin", getOrigin(), prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Electronic", getElectronic(), prec, false);
        print_utils::vector(0, "Nuclear", getNuclear(), prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::vector(0, "Total vector", getTensor(), prec, false);
        print_utils::scalar(0, "Magnitude", "(au)", length_au, prec, false);
        print_utils::scalar(0, "         ", "(Debye)", length_db, prec, false);
        mrcpp::print::separator(0, '=', 2);
    }

private:
    mrcpp::Coord<3> origin{};
    DoubleVector nuc_tensor{DoubleVector::Zero(3)};
    DoubleVector el_tensor{DoubleVector::Zero(3)};
};
// clang-format on

} // namespace mrchem

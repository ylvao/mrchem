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
class NMRShielding final {
public:
    NMRShielding(int k, const Nucleus &n) : K(k), nuc(n) {}

    int getK() const { return this->K; }
    mrcpp::Coord<3> &getOrigin() { return this->origin; }
    const mrcpp::Coord<3> &getOrigin() const { return this->origin; }
    const Nucleus &getNucleus() const { return this->nuc; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    void print() const {
        auto iso_ppm_d = getDiamagnetic().trace() / 3.0;
        auto iso_ppm_p = getParamagnetic().trace() / 3.0;
        auto iso_ppm_t = iso_ppm_d + iso_ppm_p;

        std::stringstream o_nucleus;
        o_nucleus << " Nucleus ";
        o_nucleus << std::setw(24) << getK();
        o_nucleus << std::setw(13) << getNucleus().getElement().getSymbol();

        auto prec = mrcpp::Printer::getPrecision();
        mrcpp::print::header(0, "NMR shielding");
        println(0, o_nucleus.str());
        print_utils::coord(0, "r_K", getNucleus().getCoord(), prec, false);
        print_utils::coord(0, "r_O", getOrigin(), prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Total tensor", getTensor(), prec, false);
        print_utils::scalar(0, "Iso. average", "(ppm)", iso_ppm_t, prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Diamagnetic ", getDiamagnetic(), prec, false);
        print_utils::scalar(0, "Iso. average", "(ppm)", iso_ppm_d, prec, false);
        mrcpp::print::separator(0, '-');
        print_utils::matrix(0, "Paramagnetic", getParamagnetic(), prec, false);
        print_utils::scalar(0, "Iso. average", "(ppm)", iso_ppm_d, prec, false);
        mrcpp::print::separator(0, '=', 2);
    }

private:
    const int K;
    const Nucleus nuc;
    mrcpp::Coord<3> origin{};
    DoubleMatrix dia_tensor{DoubleMatrix::Zero(3,3)};
    DoubleMatrix para_tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} // namespace mrchem

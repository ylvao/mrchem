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

    friend std::ostream& operator<<(std::ostream &o, const NMRShielding &nmr) {
        auto prec = mrcpp::Printer::getPrecision();
        auto isoDSppm = nmr.getDiamagnetic().trace() / 3.0;
        auto isoPSppm = nmr.getParamagnetic().trace() / 3.0;
        auto isoTSppm = isoDSppm + isoPSppm;

        std::stringstream o_tot_ppm, o_dia_ppm, o_para_ppm;
        o_tot_ppm << std::setw(27) << std::setprecision(prec) << std::fixed << isoTSppm;
        o_dia_ppm << std::setw(27) << std::setprecision(prec) << std::fixed << isoDSppm;
        o_para_ppm << std::setw(27) << std::setprecision(prec) << std::fixed << isoPSppm;

        std::string dia_str_0 = math_utils::vector_to_string(prec, 14, nmr.getDiamagnetic().row(0));
        std::string dia_str_1 = math_utils::vector_to_string(prec, 14, nmr.getDiamagnetic().row(1));
        std::string dia_str_2 = math_utils::vector_to_string(prec, 14, nmr.getDiamagnetic().row(2));

        std::string para_str_0 = math_utils::vector_to_string(prec, 14, nmr.getParamagnetic().row(0));
        std::string para_str_1 = math_utils::vector_to_string(prec, 14, nmr.getParamagnetic().row(1));
        std::string para_str_2 = math_utils::vector_to_string(prec, 14, nmr.getParamagnetic().row(2));

        std::string tot_str_0 = math_utils::vector_to_string(prec, 14, nmr.getTensor().row(0));
        std::string tot_str_1 = math_utils::vector_to_string(prec, 14, nmr.getTensor().row(1));
        std::string tot_str_2 = math_utils::vector_to_string(prec, 14, nmr.getTensor().row(2));

        std::stringstream o_nucleus;
        o_nucleus << std::setw(14) << nmr.getK() << std::setw(14) << nmr.getNucleus().getElement().getSymbol();

        const auto &coord = nmr.getNucleus().getCoord();
        std::string coord_str = math_utils::coord_to_string(prec, 14, coord);
        std::string origin_str = math_utils::coord_to_string(prec, 14, nmr.getOrigin());

        o << "============================================================" << std::endl;
        o << "                       NMR shielding                        " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "        Nucleus :" << o_nucleus.str()                         << std::endl;
        o << "            r_K :" << coord_str                               << std::endl;
        o << "            r_O :" << origin_str                              << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Total tensor :" << tot_str_0                               << std::endl;
        o << "                :" << tot_str_1                               << std::endl;
        o << "                :" << tot_str_2                               << std::endl;
        o << "   Iso. average :         (ppm) " << o_tot_ppm.str()          << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "    Diamagnetic :" << dia_str_0                               << std::endl;
        o << "                :" << dia_str_1                               << std::endl;
        o << "                :" << dia_str_2                               << std::endl;
        o << "   Dia. average :         (ppm) " << o_dia_ppm.str()          << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Paramagnetic :" << para_str_0                              << std::endl;
        o << "                :" << para_str_1                              << std::endl;
        o << "                :" << para_str_2                              << std::endl;
        o << "  Para. average :         (ppm) " << o_para_ppm.str()         << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        return o;
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

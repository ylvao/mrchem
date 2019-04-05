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
    mrcpp::Coord<3> &getOrigin() { return this->origin; }
    const mrcpp::Coord<3> &getOrigin() const { return this->origin; }

    DoubleMatrix getTensor() const { return getDiamagnetic() + getParamagnetic(); }
    DoubleMatrix &getDiamagnetic() { return this->dia_tensor; }
    DoubleMatrix &getParamagnetic() { return this->para_tensor; }
    const DoubleMatrix &getDiamagnetic() const { return this->dia_tensor; }
    const DoubleMatrix &getParamagnetic() const { return this->para_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const Magnetizability &mag) {
        auto prec = mrcpp::Printer::getPrecision();
        auto isoDMau = mag.getDiamagnetic().trace() / 3.0;
        auto isoPMau = mag.getParamagnetic().trace() / 3.0;
        auto isoTMau = isoDMau + isoPMau;

        std::string origin_str = math_utils::coord_to_string(prec, 14, mag.getOrigin());

        std::stringstream o_tot_au, o_dia_au, o_para_au;
        o_tot_au << std::setw(27) << std::setprecision(prec) << std::fixed << isoTMau;
        o_dia_au << std::setw(27) << std::setprecision(prec) << std::fixed << isoDMau;
        o_para_au << std::setw(27) << std::setprecision(prec) << std::fixed << isoPMau;

        // SI units (J/T^2 10^{-30})
        std::stringstream o_tot_si, o_dia_si, o_para_si;
        o_tot_si << std::setw(27) << std::setprecision(prec) << std::fixed << isoTMau * PHYSCONST::JT_m2;
        o_dia_si << std::setw(27) << std::setprecision(prec) << std::fixed << isoDMau * PHYSCONST::JT_m2;
        o_para_si << std::setw(27) << std::setprecision(prec) << std::fixed << isoPMau * PHYSCONST::JT_m2;

        std::string dia_str_0 = math_utils::vector_to_string(prec, 14, mag.getDiamagnetic().row(0));
        std::string dia_str_1 = math_utils::vector_to_string(prec, 14, mag.getDiamagnetic().row(1));
        std::string dia_str_2 = math_utils::vector_to_string(prec, 14, mag.getDiamagnetic().row(2));

        std::string para_str_0 = math_utils::vector_to_string(prec, 14, mag.getParamagnetic().row(0));
        std::string para_str_1 = math_utils::vector_to_string(prec, 14, mag.getParamagnetic().row(1));
        std::string para_str_2 = math_utils::vector_to_string(prec, 14, mag.getParamagnetic().row(2));

        std::string tot_str_0 = math_utils::vector_to_string(prec, 14, mag.getTensor().row(0));
        std::string tot_str_1 = math_utils::vector_to_string(prec, 14, mag.getTensor().row(1));
        std::string tot_str_2 = math_utils::vector_to_string(prec, 14, mag.getTensor().row(2));

        o << "============================================================" << std::endl;
        o << "                      Magnetizability                       " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "            r_O :" << origin_str                              << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Total tensor :" << tot_str_0                               << std::endl;
        o << "                :" << tot_str_1                               << std::endl;
        o << "                :" << tot_str_2                               << std::endl;
        o << "   Iso. average :          (au) " << o_tot_au.str()           << std::endl;
        o << "                :          (SI) " << o_tot_si.str()           << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "    Diamagnetic :" << dia_str_0                               << std::endl;
        o << "                :" << dia_str_1                               << std::endl;
        o << "                :" << dia_str_2                               << std::endl;
        o << "   Dia. average :          (au) " << o_dia_au.str()           << std::endl;
        o << "                :          (SI) " << o_dia_si.str()           << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "   Paramagnetic :" << para_str_0                              << std::endl;
        o << "                :" << para_str_1                              << std::endl;
        o << "                :" << para_str_2                              << std::endl;
        o << "  Para. average :          (au) " << o_para_au.str()          << std::endl;
        o << "                :          (SI) " << o_para_si.str()          << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;

        return o;
    }

private:
    mrcpp::Coord<3> origin{};
    DoubleMatrix dia_tensor{DoubleMatrix::Zero(3,3)};
    DoubleMatrix para_tensor{DoubleMatrix::Zero(3,3)};
};
// clang-format on

} // namespace mrchem

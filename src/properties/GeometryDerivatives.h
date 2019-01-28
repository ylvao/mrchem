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
class GeometryDerivatives final {
public:
    explicit GeometryDerivatives(int k)
        : nuc_tensor{DoubleMatrix::Zero(k, 3)}
        , el_tensor{DoubleMatrix::Zero(k, 3)} {}

    DoubleMatrix getTensor() const { return getNuclear() + getElectronic(); }
    DoubleMatrix &getNuclear() { return this->nuc_tensor; }
    DoubleMatrix &getElectronic() { return this->el_tensor; }
    const DoubleMatrix &getNuclear() const { return this->nuc_tensor; }
    const DoubleMatrix &getElectronic() const { return this->el_tensor; }

    friend std::ostream& operator<<(std::ostream &o, const GeometryDerivatives &gd) {
        auto length_au = gd.getTensor().norm();

        auto oldPrec = mrcpp::Printer::setPrecision(10);
        Eigen::IOFormat clean_format(10, 0, ", ", "\n", " [", "] ");
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                   Geometry derivatives                     " << std::endl;
        o << "------------------------------------------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o << " Length of vector:   (au)    " << std::setw(30) << length_au  << std::endl;
        o << "                                                            " << std::endl;
        o << "-------------------------- Total ---------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<             gd.getTensor().format(clean_format)                << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------- Nuclear --------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<             gd.getNuclear().format(clean_format)               << std::endl;
        o << "                                                            " << std::endl;
        o << "------------------------ Electronic ------------------------" << std::endl;
        o << "                                                            " << std::endl;
        o <<             gd.getElectronic().format(clean_format)            << std::endl;
        o << "                                                            " << std::endl;
        o << "============================================================" << std::endl;
        o << "                                                            " << std::endl;
        mrcpp::Printer::setPrecision(oldPrec);

        return o;
    }

private:
    DoubleMatrix nuc_tensor;
    DoubleMatrix el_tensor;
};
// clang-format on

} //namespace mrchem

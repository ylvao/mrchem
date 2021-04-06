/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/DistanceOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/PositionOperator.h"

namespace mrchem {

class X_rm3 final : public RankOneTensorOperator<3> {
public:
    X_rm3(const Nuclei &nucs, const mrcpp::Coord<3> &R_k)
            : r_m3(3.0, R_k, 1.0e-3)
            , r(nucs[0].getCoord()) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_m3 * r[0];
        h[1] = r_m3 * r[1];
        h[2] = r_m3 * r[2];
        h[0].name() = "x/r^3";
        h[1].name() = "y/r^3";
        h[2].name() = "z/r^3";
    }

private:
    DistanceOperator r_m3;
    PositionOperator r;
};

} // namespace mrchem

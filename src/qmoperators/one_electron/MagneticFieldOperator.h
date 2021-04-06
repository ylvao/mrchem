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

#include "H_B_dip.h"
#include "RankZeroTensorOperator.h"

/** @class MagneticFieldOperator
 *
 * @brief External magnetic field operator
 *
 * An external magnetic field interacts with the molecular dipole
 * moment. The operator is simply implemented as a scalar product of
 * the dipole moment operator with the magnetic field vector.
 *
 */

namespace mrchem {

class MagneticFieldOperator final : public ExternalFieldOperator {
public:
    MagneticFieldOperator(const Eigen::Vector3d &f,
                          std::shared_ptr<mrcpp::DerivativeOperator<3>> D,
                          const mrcpp::Coord<3> &o)
            : field(f)
            , dipole(D, o) {
        RankZeroTensorOperator &d_x = this->dipole[0];
        RankZeroTensorOperator &d_y = this->dipole[1];
        RankZeroTensorOperator &d_z = this->dipole[2];

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &HMF = (*this);
        HMF = f[0] * d_x + f[1] * d_y + f[2] * d_z;
        HMF.name() = "B . mu_B";
    }

private:
    Eigen::Vector3d field;
    H_B_dip dipole;
};

} // namespace mrchem

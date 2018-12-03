/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2018 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "qmfunctions/qmfunction_fwd.h"
#include "qmoperators/qmoperator_fwd.h"
#include "qmoperators/two_electron/FockOperator.h"

/** @class GroundStateHelmholtz
 *
 * @brief Container of HelmholtzOperators for a corresponding OrbtialVector
 *
 * This class assigns one HelmholtzOperator to each orbital in an OrbitalVector.
 * The operators can be re-used if several orbitals share the same (or similar)
 * energy, or if the change in energy is small relative to the previous iteration.
 */

namespace mrchem {

class GroundStateHelmholtz {
public:
    GroundStateHelmholtz(double build)
            : build_prec(build) {}

    void setup(double prec) { this->apply_prec = prec; }
    void clear() { this->apply_prec = -1.0; }

    OrbitalVector operator()(FockOperator &fock, const ComplexMatrix &F, OrbitalVector &Phi) const;

private:
    double build_prec; ///< Precision for construction of Helmholtz operators
    double apply_prec; ///< Precision for application of Helmholtz operators

    Orbital apply(mrcpp::HelmholtzOperator &H, Orbital &inp) const;
};

} // namespace mrchem

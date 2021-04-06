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

#include "MRCPP/Printer"

#include "IdentityOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** Identity operator is a deep copy */
Orbital QMIdentity::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Orbital out = inp.paramCopy();
    qmfunction::deep_copy(out, inp);
    return out;
}

/** Identity operator is a deep copy */
Orbital QMIdentity::dagger(Orbital inp) {
    return apply(inp);
}

/** Overwrite default deep copy by more efficient dot product */
ComplexDouble IdentityOperator::operator()(Orbital bra, Orbital ket) {
    if (this->I->prec() < 0.0) MSG_ERROR("Uninitialized operator");
    return orbital::dot(bra, ket);
}

/** Overwrite default deep copy by more efficient dot product */
ComplexDouble IdentityOperator::dagger(Orbital bra, Orbital ket) {
    return operator()(bra, ket);
}

/** Overwrite default deep copy by calculation of overlap matrix */
ComplexMatrix IdentityOperator::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    if (this->I->prec() < 0.0) MSG_ERROR("Uninitialized operator");

    ComplexMatrix S;
    if (&bra == &ket) {
        // In MPI it is beneficial to assume <bra| == |ket>
        S = orbital::calc_overlap_matrix(bra);
    } else {
        S = orbital::calc_overlap_matrix(bra, ket);
    }
    return S;
}

/** Overwrite default deep copy by calculation of overlap matrix */
ComplexMatrix IdentityOperator::dagger(OrbitalVector &bra, OrbitalVector &ket) {
    return operator()(bra, ket);
}

} // namespace mrchem

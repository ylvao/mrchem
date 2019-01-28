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

#include "GroundStateSolver.h"

/** @class OrbitalOptimizer
 *
 * @brief Ground state SCF optimization with kinetic energy operator
 *
 * This ground state SCF solver computes the Fock matrix explicitly by application of
 * the kinetic energy operator. This simplifies the algorithm significanly when a KAIN
 * iterative accelerator is present, or when diagonalization/localization occurs.
 * The derivative operator in the kinetic energy MIGHT affect the quadratic precision
 * in the energy, but this point seems to no longer be critical.
 */

namespace mrchem {

class Accelerator;

class OrbitalOptimizer final : public GroundStateSolver {
public:
    OrbitalOptimizer(Accelerator *k = nullptr);

    void setup(FockOperator &F, OrbitalVector &Phi, ComplexMatrix &F_mat);
    void clear();

    bool optimize();

protected:
    Accelerator *kain; ///< KAIN accelerator(pointer to external object)

    bool useKAIN() const { return (this->kain == nullptr) ? false : true; }
};

} // namespace mrchem

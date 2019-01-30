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

#include "SCF.h"
#include "properties/SCFEnergy.h"

/** @class GroundStateSolver
 *
 * @brief Abstract class for different types of ground state SCF solvers
 *
 * The ground state SCF solvers share some common features which are collected in this
 * abstract base class. This is mainly the construction of the argument that is used
 * for the Helmholtz operators.
 */

namespace mrchem {

class FockOperator;

class GroundStateSolver : public SCF {
protected:
    std::vector<SCFEnergy> energy;

    ComplexMatrix *fMat_n{nullptr};     ///< Fock matrix (pointer to external object)
    FockOperator *fOper_n{nullptr};     ///< Fock operator (pointer to external object)
    OrbitalVector *orbitals_n{nullptr}; ///< Orbtials (pointer to external object)

    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;
};

} // namespace mrchem

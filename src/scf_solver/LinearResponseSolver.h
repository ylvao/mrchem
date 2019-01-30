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

/** @class LinearResponseSolver
 *
 */

namespace mrchem {

class FockOperator;
class Accelerator;

class LinearResponseSolver final : public SCF {
public:
    LinearResponseSolver(Accelerator *k_x = nullptr, Accelerator *k_y = nullptr);

    void setupUnperturbed(double prec, FockOperator *fock, OrbitalVector *Phi, ComplexMatrix *F);
    void clearUnperturbed();

    void setup(FockOperator *fock, OrbitalVector *X);
    void setup(double omega, FockOperator *fock, OrbitalVector *X, OrbitalVector *Y);
    void clear();

    bool optimize();

protected:
    bool dynamic{false};
    double frequency{0.0};

    FockOperator *fOper_0{nullptr};
    FockOperator *fOper_1{nullptr};

    ComplexMatrix *fMat_0{nullptr};
    ComplexMatrix *fMat_x{nullptr};
    ComplexMatrix *fMat_y{nullptr};

    OrbitalVector *orbitals_0{nullptr};
    OrbitalVector *orbitals_x{nullptr};
    OrbitalVector *orbitals_y{nullptr};

    Accelerator *kain_x;
    Accelerator *kain_y;

    OrbitalVector setupHelmholtzArguments(OrbitalVector &dPhi, const ComplexMatrix &M, bool adjoint);

    void printProperty() const;
    double calcProperty();
};

} // namespace mrchem

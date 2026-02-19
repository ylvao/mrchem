
/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "chemistry/Nucleus.h"
#include "qmoperators/QMPotential.h"
#include "tensor/RankOneOperator.h"
#include "tensor/RankZeroOperator.h"
#include <string>

/** @class DFTD3Correction
 *
 * @brief Calculates D3 dispersion correction
 */

namespace mrchem {

class SCFEnergy;
class MomentumOperator;
class KineticOperator;
class ZoraKineticOperator;
class ZoraOperator;
class AZoraPotential;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;
class ReactionOperator;

class DFTD3Correction final {
public:

private:

};

} // namespace mrchem

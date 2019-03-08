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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "GroundStateSolver.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/two_electron/FockOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Computes the SCF energy by tracing the Fock operator
 *
 * Prints the current nuclear, electronic and total energies.
 */
double GroundStateSolver::calcProperty(FockOperator &F, OrbitalVector &Phi, ComplexMatrix &F_mat) {
    Printer::printHeader(0, "Calculating SCF energy");
    Timer timer;

    SCFEnergy E = F.trace(Phi, F_mat);
    this->energy.push_back(E);

    timer.stop();
    int oldPrec = Printer::setPrecision(15);
    println(0, " Nuclear energy              " << std::setw(30) << E.getNuclearEnergy());
    println(0, " Electronic energy           " << std::setw(30) << E.getElectronicEnergy());
    Printer::printSeparator(0, '-');
    println(0, " Total energy                " << std::setw(30) << E.getTotalEnergy());
    Printer::printFooter(0, timer, 2);
    Printer::setPrecision(oldPrec);

    return E.getTotalEnergy();
}

/** @brief Computes the SCF energy update from last iteration */
double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return std::abs(getUpdate(this->property, iter, false));
}

/** @brief Pretty printing of the different contributions to the SCF energy */
void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double phi_0 = scf_0.getOrbitalEnergy();
    double phi_1 = scf_1.getOrbitalEnergy();
    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();
    double N_0 = scf_0.getNuclearEnergy();
    double N_1 = scf_1.getNuclearEnergy();

    Printer::printHeader(0, "                    Energy                 Update      Done ");
    printUpdate(" Orbital    ", phi_1, phi_1 - phi_0);
    printUpdate(" Kinetic    ", T_1, T_1 - T_0);
    printUpdate(" N-E        ", V_1, V_1 - V_0);
    printUpdate(" Coulomb    ", J_1, J_1 - J_0);
    printUpdate(" Exchange   ", K_1, K_1 - K_0);
    printUpdate(" X-C        ", XC_1, XC_1 - XC_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Electronic ", E_1, E_1 - E_0);
    printUpdate(" Nuclear    ", N_1, N_1 - N_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Total      ", E_1 + N_1, (E_1 + N_1) - (E_0 + N_0));
    Printer::printSeparator(0, '=');
}

/** @brief Reset accumulated data */
void GroundStateSolver::reset() {
    SCFSolver::reset();
    this->energy.clear();
}

} // namespace mrchem

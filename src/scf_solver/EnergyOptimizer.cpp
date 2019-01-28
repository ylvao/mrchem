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

#include "EnergyOptimizer.h"
#include "HelmholtzVector.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Prepare solver for optimization
 *
 * @param fock: FockOperator defining the SCF equations
 * @param Phi: Orbitals to optimize
 * @param F: Fock matrix
 * @param fock_np1: Next iteration FockOperator
 * @param Phi_np1: Next iteration orbitals
 *
 * SCF solver will NOT take ownership of the input, so these objects must be taken
 * care of externally (do not delete until SCF goes out of scope).
 */
void EnergyOptimizer::setup(FockOperator &fock,
                            OrbitalVector &Phi,
                            ComplexMatrix &F,
                            FockOperator &fock_np1,
                            OrbitalVector &Phi_np1) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->fOper_np1 = &fock_np1;
    this->orbitals_n = &Phi;
    this->orbitals_np1 = &Phi_np1;
}

/** @brief Clear solver after optimization
 *
 * Clear pointers that was set during setup, and reset the precision parameter
 * (only the current precision orbPrec[0], not the bounaary values orbPrec[1,2]).
 * Solver can be re-used after another setup.
 */
void EnergyOptimizer::clear() {
    this->fMat_n = nullptr;
    this->fOper_n = nullptr;
    this->fOper_np1 = nullptr;
    this->orbitals_n = nullptr;
    this->orbitals_np1 = nullptr;
    resetPrecision();
}

/** @brief Run energy optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm does NOT
 * use any kinetic energy operator, but computes the Fock matrix as a potential
 * update from the previous iteration. No KAIN accelerator and no diagonalization
 * or localization within the SCF iterations. Main points of the algorithm:
 *
 * Pre SCF: diagonalize/localize orbitals
 *
 * 1) Setup Fock operator
 * 2) Compute current SCF energy
 * 3) Apply Helmholtz operator on all orbitals
 * 4) Compute orbital updates
 * 5) Compute errors and check for convergence
 * 6) Compute Fock matrix update
 * 7) Orthonormalize orbitals (Löwdin)
 *
 * Post SCF: diagonalize/localize orbitals
 */
bool EnergyOptimizer::optimize() {
    ComplexMatrix &F_mat_n = *this->fMat_n;
    OrbitalVector &Phi_n = *this->orbitals_n;
    OrbitalVector &Phi_np1 = *this->orbitals_np1;
    FockOperator &F = *this->fOper_n;
    RankZeroTensorOperator &V = F.potential();

    double orb_prec = this->orbPrec[0];
    double err_o = orbital::get_errors(Phi_n).maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    if (this->canonical) {
        orbital::diagonalize(orb_prec, Phi_n, F_mat_n);
    } else {
        ComplexMatrix U_mat = orbital::localize(orb_prec, Phi_n);
        F_mat_n = U_mat * F_mat_n * U_mat.adjoint();
    }

    int nIter = 0;
    bool converged = false;
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        // Initialize SCF cycle
        Timer timer;
        printCycleHeader(nIter);
        orb_prec = adjustPrecision(err_o);

        // Compute electronic energy
        F.setup(orb_prec);
        double E = calcProperty();
        this->property.push_back(E);

        // Apply Helmholtz operator
        HelmholtzVector H(orb_prec, F_mat_n.real().diagonal());
        ComplexMatrix L_mat = H.getLambdaMatrix().cast<ComplexDouble>();
        OrbitalVector Psi = orbital::rotate(L_mat - F_mat_n, Phi_n);
        Phi_np1 = H(V, Phi_n, Psi);
        Psi.clear();

        // Compute orbital updates
        OrbitalVector dPhi_n = orbital::add(1.0, Phi_np1, -1.0, Phi_n);

        // Compute orbital errors
        DoubleVector errors = orbital::get_norms(dPhi_n);

        err_o = errors.maxCoeff();
        err_t = errors.norm();
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);
        converged = checkConvergence(err_o, err_p);

        // Compute Fock matrix
        ComplexMatrix F_mat_np1 = F_mat_n + calcFockMatrixUpdate(orb_prec, dPhi_n, L_mat);
        dPhi_n.clear();
        F.clear();

        // Rotate orbitals
        ComplexMatrix U_mat = orbital::calc_lowdin_matrix(Phi_np1);
        Phi_n = orbital::rotate(U_mat, Phi_np1, orb_prec);
        F_mat_n = U_mat * F_mat_np1 * U_mat.adjoint();
        Phi_np1.clear();

        orbital::set_errors(Phi_n, errors);

        timer.stop();
        printOrbitals(F_mat_n.real().diagonal(), Phi_n, 0);
        printProperty();
        printCycleFooter(timer.getWallTime());

        if (converged) break;
    }

    if (this->canonical) {
        orbital::diagonalize(orb_prec / 10, Phi_n, F_mat_n);
    } else {
        ComplexMatrix U_mat = orbital::localize(orb_prec / 10, Phi_n);
        F_mat_n = U_mat * F_mat_n * U_mat.adjoint();
    }

    printConvergence(converged);
    return converged;
}

/** @brief Compute Fock matrix update
 *
 * The Fock matrix is computed as a direct update from the previous iteration.
 * This update is exact, provided that the orbital relation between the orbitals
 * at iterations n and n+1 is exactly the application of the Helmholtz operator.
 *
 * \Delta F = \Delta S_1 F + \Delta S_2 \Lambda + \Delta F_{pot}
 *
 * (\Delta S_1^n)_{ij} = <\Delta \phi_i^{n}   | \phi_j^n>
 * (\Delta S_2^n)_{ij} = <\Delta \phi_i^{n+1} | \Delta \phi_j^n>
 *
 * (\Delta F_^n{pot})_{ij} = <\phi_i^{n+1} | V^n | \Delta \phi_j^n>
 *                         + <\phi_i^{n+1} | \Delta V^n | \phi_j^n>
 *
 */
ComplexMatrix EnergyOptimizer::calcFockMatrixUpdate(double prec, OrbitalVector &dPhi_n, const ComplexMatrix &L) {
    if (this->fOper_np1 == nullptr) MSG_FATAL("Operator not initialized");

    OrbitalVector &Phi_n = *this->orbitals_n;
    OrbitalVector &Phi_np1 = *this->orbitals_np1;

    Printer::printHeader(0, "Computing Fock matrix update");

    Timer timer;
    ComplexMatrix dS_1 = orbital::calc_overlap_matrix(dPhi_n, Phi_n);
    ComplexMatrix dS_2 = orbital::calc_overlap_matrix(Phi_np1, dPhi_n);

    NuclearOperator *v_n = this->fOper_n->getNuclearOperator();
    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
    XCOperator *xc_n = this->fOper_n->getXCOperator();

    ComplexMatrix dV_n;
    { // Nuclear potential matrix is computed explicitly
        Timer timer;
        dV_n = (*v_n)(Phi_np1, dPhi_n);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Nuclear potential matrix", t, 5);
    }

    ComplexMatrix F_n;
    { // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_n(nullptr, nullptr, j_n, k_n, xc_n);
        fock_n.build();
        F_n = fock_n(Phi_np1, Phi_n);
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Fock matrix n", t, 5);
    }

    { // The n+1 Fock operator needs orthonormalized orbitals
        orbital::orthonormalize(prec, Phi_np1);
    }

    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

    println(0, "                                                            ");
    // Do not setup internal exchange, it must be applied on the fly anyway
    if (j_np1 != nullptr) j_np1->setup(prec);
    if (k_np1 != nullptr) k_np1->setup(prec);
    if (xc_np1 != nullptr) xc_np1->setup(prec);
    println(0, "                                                            ");

    ComplexMatrix F_np1;
    { // Computing potential matrix excluding nuclear part
        Timer timer;
        FockOperator fock_np1(nullptr, nullptr, j_np1, k_np1, xc_np1);
        fock_np1.build();
        ComplexMatrix F_1 = fock_np1(Phi_n, Phi_n);
        ComplexMatrix F_2 = fock_np1(Phi_n, dPhi_n);

        F_np1 = F_1 + F_2 + F_2.transpose();
        // ComplexMatrix F_3 = f_np1(*this->dPhi_n, *this->phi_n);
        // ComplexMatrix F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
        // ComplexMatrix F_np1 = F_1 + F_2 + F_3 + F_4;
        timer.stop();
        double t = timer.getWallTime();
        Printer::printDouble(0, "Fock matrix n+1", t, 5);
    }
    if (j_np1 != 0) j_np1->clear();
    if (k_np1 != 0) k_np1->clear();
    if (xc_np1 != 0) xc_np1->clear();

    // Re-computing non-orthogonal phi_np1
    Phi_np1 = orbital::add(1.0, Phi_n, 1.0, dPhi_n);

    ComplexMatrix dF_1 = dS_1 * (*this->fMat_n);
    ComplexMatrix dF_2 = dS_2 * L;
    ComplexMatrix dF_3 = F_np1 - F_n;

    // Adding up the pieces
    ComplexMatrix dF_n = dV_n + dF_1 + dF_2 + dF_3;

    // Symmetrizing
    ComplexMatrix sym = dF_n + dF_n.transpose();
    dF_n = 0.5 * sym;

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return dF_n;
}

} // namespace mrchem

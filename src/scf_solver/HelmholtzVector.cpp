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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "HelmholtzVector.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief HelmholtzVector constructor
 *
 * This will set the build precision of the Helmholtz operators and the vector
 * of lambda parameters that will be used in the subsequent application. No
 * operators are constructed at this point, they are produced on-the-fly in
 * the application.
 */
HelmholtzVector::HelmholtzVector(double pr, const DoubleVector &l)
        : prec(pr) {
    this->lambda = l;
    for (int i = 0; i < this->lambda.size(); i++) {
        if (this->lambda(i) > 0.0) this->lambda(i) = -0.5;
    }
}

/** @brief Apply Helmholtz operator component wise on OrbitalVector
 *
 * This will construct a separate Helmholtz operator for each of the entries
 * in the OrbitalVector based on the corresponding lambda_i parameter in the
 * HelmholtzVector. Computes output as: out_i = H_i[phi_i]
 *
 * MPI: Output vector gets the same MPI distribution as input vector. Only
 *      local orbitals are computed.
 */
OrbitalVector HelmholtzVector::operator()(OrbitalVector &Phi) const {
    mrcpp::print::header(0, "Applying Helmholtz operators");

    Timer t_tot;
    int pprec = Printer::getPrecision();
    OrbitalVector out = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;

        Timer t_i;
        out[i] = apply(i, Phi[i]);
        out[i].rescale(-1.0 / (2.0 * MATHCONST::pi));

        std::stringstream o_txt;
        o_txt << std::setw(4) << i;
        o_txt << std::setw(19) << std::setprecision(pprec) << std::scientific << out[i].norm();
        print_utils::qmfunction(0, o_txt.str(), out[i], t_i);
    }
    mrcpp::print::footer(0, t_tot, 2);
    return out;
}

/** @brief Prepare the matrix term in the Helmholtz argument
 *
 * This will subtract the lambda parameters in the Helmholtz operators from the
 * diagonal of the Fock matrix and rotate the orbitals with the resulting matrix
 * psi_i = \sum_j(L_ij - F_ij)*phi_j
 *
 * MPI: Output vector gets the same MPI distribution as input vector. Only
 *      local orbitals are computed.
 */
OrbitalVector HelmholtzVector::rotate(const ComplexMatrix &F_mat, OrbitalVector &Phi) const {
    Timer t_tot;
    mrcpp::print::header(0, "Rotating Helmholtz argument");
    ComplexMatrix L_mat = getLambdaMatrix();

    Timer rot_t;
    OrbitalVector Psi = orbital::rotate(L_mat - F_mat, Phi);
    mrcpp::print::time(0, "Rotating orbitals", rot_t);

    mrcpp::print::footer(0, t_tot, 2);
    return Psi;
}

/** @brief Apply Helmholtz operator component wise on OrbitalVector
 *
 * This will construct a separate Helmholtz operator for each of the entries
 * in the OrbitalVector based on the corresponding lambda_i parameter in the
 * HelmholtzVector. Computes output as: out_i = H_i[V*phi_i + psi_i]
 *
 * Specialized version with smaller memory footprint since the full vector V*Phi
 * is never stored, but computed on the fly.
 *
 * MPI: Output vector gets the same MPI distribution as input vector. Only
 *      local orbitals are computed.
 */
OrbitalVector HelmholtzVector::apply(RankZeroTensorOperator &V, OrbitalVector &Phi, OrbitalVector &Psi) const {
    mrcpp::print::header(0, "Applying Helmholtz operators");
    if (Phi.size() != Psi.size()) MSG_ABORT("OrbitalVector size mismatch");

    Timer t_tot;
    int pprec = Printer::getPrecision();
    OrbitalVector out = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;

        Timer t_i;
        Orbital Vphi_i = V(Phi[i]);
        Vphi_i.add(1.0, Psi[i]);
        Vphi_i.rescale(-1.0 / (2.0 * MATHCONST::pi));
        out[i] = apply(i, Vphi_i);

        std::stringstream o_txt;
        o_txt << std::setw(4) << i;
        o_txt << std::setw(19) << std::setprecision(pprec) << std::scientific << out[i].norm();
        print_utils::qmfunction(0, o_txt.str(), out[i], t_i);
    }
    mrcpp::print::footer(0, t_tot, 2);
    return out;
}

/** @brief Apply Helmholtz operator on individual Orbital
 *
 * This will construct a Helmholtz operator with the i-th component of the
 * lambda vector and apply it to the input orbital.
 *
 * Computes output as: out_i = H_i[phi_i]
 */
Orbital HelmholtzVector::apply(int i, Orbital &phi) const {
    ComplexDouble mu_i = std::sqrt(-2.0 * this->lambda(i));
    if (std::abs(mu_i.imag()) > mrcpp::MachineZero) MSG_ABORT("Mu cannot be complex");
    mrcpp::HelmholtzOperator H(*MRA, mu_i.real(), this->prec);

    Orbital out = phi.paramCopy();
    if (phi.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(this->prec, out.real(), H, phi.real());
    }
    if (phi.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(this->prec, out.imag(), H, phi.imag());
        if (phi.conjugate()) out.imag().rescale(-1.0);
    }
    return out;
}
} // namespace mrchem

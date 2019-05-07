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

#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "core.h"
#include "gto.h"
#include "parallel.h"
#include "sad.h"
#include "utils/math_utils.h"

#include "chemistry/Molecule.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {

ComplexMatrix diagonalize_fock(KineticOperator &T, RankZeroTensorOperator &V, OrbitalVector &Phi, int spin);
OrbitalVector rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin);
void project_atomic_densities(double prec, Density &rho_tot, const Molecule &mol);

} // namespace sad
} // namespace initial_guess

OrbitalVector initial_guess::sad::setup(double prec, const Molecule &mol, bool restricted, int zeta) {
    std::string restr_str = "False";
    if (restricted) restr_str = "True";
    mrcpp::print::separator(0, '-');
    println(0, " Method            : Diagonalize Hamiltonian matrix");
    println(0, " Precision         : " << prec);
    println(0, " Restricted        : " << restr_str);
    println(0, " Hamiltonian       : Superposition of Atomic Densities (SAD)");
    println(0, " Functional        : LDA (SVWN5)");
    println(0, " AO basis          : Hydrogenic orbitals");
    println(0, " Zeta quality      : " << zeta);
    mrcpp::print::separator(0, '-', 2);

    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity(); // multiplicity
    int Ne = mol.getNElectrons();     // total electrons
    int Nd = Ne - (mult - 1);         // doubly occupied electrons
    if (Nd % 2 != 0) MSG_ABORT("Invalid multiplicity");
    int Na = Nd / 2 + (mult - 1); // alpha orbitals
    int Nb = Nd / 2;              // beta orbitals

    // Make Fock operator contributions
    auto P_p = std::make_shared<mrcpp::PoissonOperator>(*MRA, prec);
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    auto xcfun_p = std::make_shared<mrdft::XCFunctional>(*MRA, not(restricted));
    xcfun_p->setFunctional("SLATERX");
    xcfun_p->setFunctional("VWN5C");
    xcfun_p->evalSetup(1);
    KineticOperator T(D_p);
    NuclearOperator V_nuc(mol.getNuclei(), prec);
    CoulombOperator J(P_p);
    XCOperator XC(xcfun_p);
    RankZeroTensorOperator V = V_nuc + J + XC;

    // Compute Coulomb density
    Density &rho_j = J.getDensity();
    initial_guess::sad::project_atomic_densities(prec, rho_j, mol);

    // Compute XC density
    if (restricted) {
        mrcpp::FunctionTree<3> &rho_xc = XC.getDensity(DENSITY::Total);
        mrcpp::copy_grid(rho_xc, rho_j.real());
        mrcpp::copy_func(rho_xc, rho_j.real());
    } else {
        mrcpp::FunctionTree<3> &rho_a = XC.getDensity(DENSITY::Alpha);
        mrcpp::FunctionTree<3> &rho_b = XC.getDensity(DENSITY::Beta);
        mrcpp::add(prec, rho_a, 1.0, rho_j.real(), -1.0 * Nb / Ne, rho_j.real());
        mrcpp::add(prec, rho_b, 1.0, rho_j.real(), -1.0 * Na / Ne, rho_j.real());

        // Extend to union grid
        int nNodes = 1;
        while (nNodes > 0) {
            int nAlpha = mrcpp::refine_grid(rho_a, rho_b);
            int nBeta = mrcpp::refine_grid(rho_b, rho_a);
            nNodes = nAlpha + nBeta;
        }
    }

    // Project AO basis of hydrogen functions
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);

    mrcpp::print::header(0, "Setting up Fock operator");
    Timer t_fock;
    T.setup(prec);
    V.setup(prec);
    mrcpp::print::footer(0, t_fock, 2);

    // Compute Fock matrix
    mrcpp::print::header(0, "Diagonalize Fock matrix");
    Timer t_diag;
    OrbitalVector Psi;
    if (restricted) {
        if (mult != 1) MSG_ABORT("Restricted open-shell not available");
        int Np = Nd / 2; // paired orbitals
        ComplexMatrix U = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Paired);
        Psi = initial_guess::sad::rotate_orbitals(prec, U, Phi, Np, SPIN::Paired);
    } else {
        int Na = Nd / 2 + (mult - 1); // alpha orbitals
        int Nb = Nd / 2;              // beta orbitals

        ComplexMatrix U_a = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Alpha);
        OrbitalVector Psi_a = initial_guess::sad::rotate_orbitals(prec, U_a, Phi, Na, SPIN::Alpha);

        ComplexMatrix U_b = initial_guess::sad::diagonalize_fock(T, V, Phi, SPIN::Beta);
        OrbitalVector Psi_b = initial_guess::sad::rotate_orbitals(prec, U_b, Phi, Nb, SPIN::Beta);

        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    T.clear();
    V.clear();
    mrcpp::print::footer(0, t_diag, 2);

    return Psi;
}

ComplexMatrix initial_guess::sad::diagonalize_fock(KineticOperator &T,
                                                   RankZeroTensorOperator &V,
                                                   OrbitalVector &Phi,
                                                   int spin) {
    Timer t1;
    for (auto &i : Phi) i.setSpin(spin);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    ComplexMatrix f_tilde = T(Phi, Phi) + V(Phi, Phi);
    ComplexMatrix f = S_m12.adjoint() * f_tilde * S_m12;
    mrcpp::print::time(0, "Compute Fock matrix", t1);

    Timer t2;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(f.cols());
    es.compute(f);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    mrcpp::print::time(0, "Diagonalize Fock matrix", t2);

    return U;
}

OrbitalVector initial_guess::sad::rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin) {
    Timer timer;
    OrbitalVector Psi;
    for (int i = 0; i < N; i++) Psi.push_back(Orbital(spin));
    mpi::distribute(Psi);

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < Psi.size(); i++) {
            if (not mpi::my_orb(Psi[i])) continue;
            QMFunctionVector func_vec;
            ComplexVector coef_vec(iter.get_size());
            for (int j = 0; j < iter.get_size(); j++) {
                int idx_j = iter.idx(j);
                Orbital &recv_j = iter.orbital(j);
                coef_vec[j] = U(i, idx_j);
                func_vec.push_back(recv_j);
            }
            Orbital tmp_i = Psi[i].paramCopy();
            qmfunction::linear_combination(tmp_i, coef_vec, func_vec, prec);
            Psi[i].add(1.0, tmp_i); // In place addition
            Psi[i].crop(prec);
        }
    }
    mrcpp::print::time(0, "Rotate orbitals", timer);
    return Psi;
}

void initial_guess::sad::project_atomic_densities(double prec, Density &rho_tot, const Molecule &mol) {
    Timer timer;
    mrcpp::print::header(0, "Projecting GTO density");
    println(0, "    N    Atom          Nuclear charge          Rho_K");
    mrcpp::print::separator(0, '-');

    auto print_prec = Printer::getPrecision();
    auto crop_prec = (mpi::numerically_exact) ? -1.0 : prec;

    std::string sad_path = SAD_BASIS_DIR;

    Density rho_loc(false);
    rho_loc.alloc(NUMBER::Real);
    rho_loc.real().setZero();

    const Nuclei &nucs = mol.getNuclei();
    for (int k = 0; k < nucs.size(); k++) {
        if (mpi::orb_rank != k % mpi::orb_size) continue;

        const std::string &sym = nucs[k].getElement().getSymbol();
        std::stringstream o_bas, o_dens;
        o_bas << sad_path << sym << ".bas";
        o_dens << sad_path << sym << ".dens";

        Density rho_k = initial_guess::gto::project_density(prec, nucs[k], o_bas.str(), o_dens.str());

        std::stringstream o_nuc, o_rho;
        o_nuc << std::setprecision(2 * print_prec) << std::fixed << nucs[k].getCharge();
        o_rho << std::setprecision(2 * print_prec) << std::fixed << rho_k.integrate().real();

        printout(0, std::setw(5) << k);
        printout(0, std::setw(7) << sym);
        printout(0, std::setw(25) << o_nuc.str());
        printout(0, std::setw(22) << o_rho.str() << std::endl);

        rho_loc.add(1.0, rho_k);
        rho_loc.crop(crop_prec);
    }

    density::allreduce_density(prec, rho_tot, rho_loc);
    mrcpp::print::footer(0, timer, 2);
}

} // namespace mrchem

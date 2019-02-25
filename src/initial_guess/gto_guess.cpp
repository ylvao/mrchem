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
#include "getkw/Getkw.hpp"

#include "mrchem.h"
#include "mrenv.h"
#include "parallel.h"

#include "chemistry/Molecule.h"
#include "gto.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

Getkw mrchem::input;
nlohmann::json mrchem::json_input;

mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

using namespace mrcpp;
using namespace mrchem;

/** @file gto_guess.cpp
 *
 * Standalone executable (gto-guess) for reading a GTO initial guess and
 * writing the resulting MW orbitals to disk.
 *
 * Requires the following input files (file names can be changed in input):
 * @mrchem.inp: regular input file, parsed through getkw (./mrchem -D)
 * initial_guess/mrchem.bas: basis set file (LSDalton format)
 * initial_guess/mrchem.moa: MO matrix (alpha orbitals or closed-shell)
 * initial_guess/mrchem.mob: MO matrix (beta orbitals)
 *
 * Produces the following output files (file names can be changed in input):
 * orbitals/phi_0.meta: orbital meta data
 * orbitals/phi_0_re.tree: MW representation of real part
 * orbitals/phi_1.meta: orbital meta data
 * orbitals/phi_1_re.tree: MW representation of real part
 */

int main(int argc, char **argv) {
    mpi::initialize(argc, argv);
    mrenv::initialize(argc, argv);
    Timer timer;

    // Reading input
    double prec = input.get<double>("rel_prec");
    bool wf_restricted = input.get<bool>("wavefunction.restricted");
    int mol_charge = input.get<int>("molecule.charge");
    int mol_multiplicity = input.get<int>("molecule.multiplicity");
    std::vector<std::string> mol_coords = input.getData("molecule.coords");
    std::string scf_guess = input.get<std::string>("scf.initial_guess");
    std::string orb_file = input.get<std::string>("files.start_orbitals");
    std::string bas_file = input.get<std::string>("files.basis_set");
    std::string moa_file = input.get<std::string>("files.mo_mat_a");
    std::string mob_file = input.get<std::string>("files.mo_mat_b");

    if (scf_guess != "gto") MSG_FATAL("Invalid initial guess");

    // Setting up molecule
    Molecule mol(mol_coords, mol_charge, mol_multiplicity);
    mol.printGeometry();

    // Setting up orbitals
    OrbitalVector Phi;
    if (wf_restricted) {
        Phi = initial_guess::gto::setup(prec, mol, bas_file, moa_file);
    } else {
        Phi = initial_guess::gto::setup(prec, mol, bas_file, moa_file, mob_file);
    }
    orbital::save_orbitals(Phi, orb_file);

    timer.stop();
    mrenv::finalize(timer.getWallTime());
    mpi::finalize();

    return 0;
}

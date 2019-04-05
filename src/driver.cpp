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

#include "driver.h"

#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"

#include "initial_guess/core.h"
#include "initial_guess/gto.h"
#include "initial_guess/sad.h"

#include "utils/math_utils.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/two_electron/CoulombOperator.h"
#include "qmoperators/two_electron/ExchangeOperator.h"
#include "qmoperators/two_electron/FockOperator.h"
#include "qmoperators/two_electron/XCOperator.h"

#include "qmoperators/one_electron/H_BB_dia.h"
#include "qmoperators/one_electron/H_BM_dia.h"
#include "qmoperators/one_electron/H_B_dip.h"
#include "qmoperators/one_electron/H_B_spin.h"
#include "qmoperators/one_electron/H_E_dip.h"
#include "qmoperators/one_electron/H_M_fc.h"
#include "qmoperators/one_electron/H_M_pso.h"
#include "qmoperators/one_electron/NuclearGradientOperator.h"

#include "scf_solver/EnergyOptimizer.h"
#include "scf_solver/KAIN.h"
#include "scf_solver/LinearResponseSolver.h"
#include "scf_solver/OrbitalOptimizer.h"

#include "mrdft/XCFunctional.h"

using mrcpp::ABGVOperator;
using mrcpp::Coord;
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

using mrdft::XCFunctional;

using nlohmann::json;

using DerivativeOperator = mrcpp::DerivativeOperator<3>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

using PoissonOperator = mrcpp::PoissonOperator;
using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;

extern mrcpp::MultiResolutionAnalysis<3> *mrchem::MRA;

namespace mrchem {

namespace driver {
void build_fock_operator(const json &input, Molecule &mol, FockOperator &F, int order);

void calc_scf_properties(const json &input, Molecule &mol);
void calc_rsp_properties(const json &input, Molecule &mol, int dir, double omega);

RankOneTensorOperator<3> get_perturbation(const json &input);
DerivativeOperator_p get_derivative(const std::string &name);
} // namespace driver

/** @brief Initialize a molecule from input
 *
 * This function expects the "molecule" subsection of the input.
 */
void driver::init_molecule(const json &json_mol, Molecule &mol) {
    println(0, "                                                            ");
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***               Initializing Molecule                  ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, "                                                            ");
    println(0, "                                                            ");

    auto charge = json_mol["charge"].get<int>();
    auto multiplicity = json_mol["multiplicity"].get<int>();

    mol.setCharge(charge);
    mol.setMultiplicity(multiplicity);

    Nuclei &nuclei = mol.getNuclei();
    for (const auto &coord : json_mol["coords"].get<json>()) {
        auto atom = coord["atom"].get<std::string>();
        auto xyz = coord["xyz"].get<std::array<double, 3>>();
        nuclei.push_back(atom, xyz);
    }

    mol.printGeometry();
}

/** @brief Run initial guess calculation for the orbitals
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule, based on the chosen initial guess method.
 * The orbital vector is initialized with the appropriate particle
 * number and spin. The Fock matrix is initialized to the zero matrix
 * of the appropriate size.
 *
 * This function expects the "initial_guess" subsection of the input.
 */
bool driver::run_guess(const json &json_guess, Molecule &mol) {
    Printer::printHeader(0, "Initial guess input");
    println(0, json_guess.dump(2));
    Printer::printSeparator(0, '=', 2);

    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    auto method = json_guess["method"].get<std::string>();
    if (method == "mw") {
        auto start_orbs = json_guess["start_orbitals"].get<std::string>();
        Phi = orbital::load_orbitals(start_orbs);
    } else if (method == "core") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto zeta = json_guess["zeta"].get<int>();
        Phi = initial_guess::core::setup(guess_prec, mol, restricted, zeta);
    } else if (method == "sad") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto zeta = json_guess["zeta"].get<int>();
        Phi = initial_guess::sad::setup(guess_prec, mol, restricted, zeta);
    } else if (method == "gto") {
        auto guess_prec = json_guess["guess_prec"].get<double>();
        auto restricted = json_guess["restricted"].get<bool>();
        auto file_basis = json_guess["file_basis"].get<std::string>();
        auto file_moa = json_guess["file_moa"].get<std::string>();
        auto file_mob = json_guess["file_mob"].get<std::string>();
        if (restricted) {
            Phi = initial_guess::gto::setup(guess_prec, mol, file_basis, file_moa);
        } else {
            Phi = initial_guess::gto::setup(guess_prec, mol, file_basis, file_moa, file_mob);
        }
    } else {
        MSG_ERROR("Invalid initial guess");
        return false;
    }

    F_mat = ComplexMatrix::Zero(Phi.size(), Phi.size());

    auto write_orbs = json_guess["write_orbitals"].get<bool>();
    auto final_orbs = json_guess["final_orbitals"].get<std::string>();
    if (write_orbs) orbital::save_orbitals(Phi, final_orbs);

    return true;
}

/** @brief Run ground-state SCF calculation
 *
 * This function will update the ground state orbitals and the Fock
 * matrix of the molecule based on the chosen electronic structure
 * method. The resulting orbitals will be either diagonalized or
 * localized at exit. Returns true if the calculation converges.
 *
 * After convergence the requested ground-state properties are computed.
 *
 * This function expects the "scf_calculation" subsection of the input.
 */
bool driver::run_scf(const json &json_scf, Molecule &mol) {
    Printer::printHeader(0, "SCF input");
    println(0, json_scf.dump(2));
    Printer::printSeparator(0, '=', 2);

    auto success = true;
    auto scf_prec = json_scf["scf_prec"].get<double>();
    auto localize = json_scf["localize"].get<bool>();

    const auto &json_fock = json_scf["fock_operator"].get<json>();
    FockOperator F;
    driver::build_fock_operator(json_fock, mol, F, 0);

    ///////////////////////////////////////////////////////////
    /////////////////   Running SCF Solver  ///////////////////
    ///////////////////////////////////////////////////////////

    auto &Phi = mol.getOrbitals();
    auto &F_mat = mol.getFockMatrix();

    // Run OrbitalOptimizer if present in input JSON
    auto orbital_solver_it = json_scf.find("orbital_solver");
    if (orbital_solver_it != json_scf.end()) {
        const auto &json_solver = *orbital_solver_it;
        auto kain = json_solver["kain"].get<int>();
        auto max_iter = json_solver["max_iter"].get<int>();
        auto rotation = json_solver["rotation"].get<int>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        OrbitalOptimizer solver;
        solver.setHistory(kain);
        solver.setMaxIterations(max_iter);
        solver.setRotation(rotation);
        solver.setCanonical(not(localize));
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);
        success = solver.optimize(F, Phi, F_mat);
    } else {
        F.setup(scf_prec);
        F_mat = F(Phi, Phi);
        F.clear();

        if (localize) {
            orbital::localize(scf_prec, Phi, F_mat);
        } else {
            orbital::diagonalize(scf_prec, Phi, F_mat);
        }
    }

    // Run EnergyOptimizer if present in input JSON
    auto energy_solver_it = json_scf.find("energy_solver");
    if (energy_solver_it != json_scf.end()) {
        const auto &json_solver = *energy_solver_it;
        auto max_iter = json_solver["max_iter"].get<int>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        EnergyOptimizer solver;
        solver.setMaxIterations(max_iter);
        solver.setRotation(1);
        solver.setCanonical(not(localize));
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);
        success = solver.optimize(F, Phi, F_mat);
    }

    ///////////////////////////////////////////////////////////
    ///////////////   Computing Properties   //////////////////
    ///////////////////////////////////////////////////////////

    F.setup(scf_prec);
    Printer::printHeader(0, "Calculating SCF energy");
    Timer timer;
    SCFEnergy &energy = mol.getSCFEnergy();
    energy = F.trace(Phi, F_mat);
    timer.stop();
    Printer::printFooter(0, timer, 2);
    F.clear();

    if (success) {
        const auto &json_prop = json_scf["properties"].get<json>();
        driver::calc_scf_properties(json_prop, mol);

        auto final_orbs = json_scf["final_orbitals"].get<std::string>();
        auto write_orbs = json_scf["write_orbitals"].get<bool>();
        if (write_orbs) orbital::save_orbitals(Phi, final_orbs);
    }

    return success;
}

/** @brief Run linear response SCF calculation
 *
 * This function will update the perturbed orbitals of the molecule
 * based on the chosen electronic structure method and perturbation
 * operator. Each response calculation corresponds to one particular
 * perturbation operator (could be a vector operator with several
 * components). Returns true if the calculation converges.
 *
 * After convergence the requested linear response properties are computed.
 *
 * This function expects a single subsection entry in the "rsp_calculations"
 * vector of the input.
 */
bool driver::run_rsp(const json &json_rsp, Molecule &mol) {
    Printer::printHeader(0, "Response input");
    println(0, json_rsp.dump(2));
    Printer::printSeparator(0, '=', 2);

    auto success = true;
    auto rsp_prec = json_rsp["rsp_prec"].get<double>();
    auto dynamic = json_rsp["dynamic"].get<bool>();
    auto localize = json_rsp["localize"].get<bool>();

    mol.initPerturbedOrbitals(dynamic);

    // Setup Fock operators
    const auto &json_fock_0 = json_rsp["fock_operator_0"].get<json>();
    const auto &json_fock_1 = json_rsp["fock_operator_1"].get<json>();
    FockOperator F_0, F_1;
    driver::build_fock_operator(json_fock_0, mol, F_0, 0);
    driver::build_fock_operator(json_fock_1, mol, F_1, 1);

    F_1.getXCOperator()->setupDensity(rsp_prec);
    F_1.getXCOperator()->setupPotential(rsp_prec);

    auto &F_mat = mol.getFockMatrix();
    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    // Setup perturbation operator
    const auto &json_pert = json_rsp["perturbation"].get<json>();
    auto h_1 = driver::get_perturbation(json_pert);

    ///////////////////////////////////////////////////////////
    /////////////////   Running RSP Solver  ///////////////////
    ///////////////////////////////////////////////////////////

    auto rsp_solver_it = json_rsp.find("rsp_solver");
    if (rsp_solver_it != json_rsp.end()) {
        const auto &json_solver = json_rsp["rsp_solver"].get<json>();
        auto kain = json_solver["kain"].get<int>();
        auto omega = json_solver["frequency"].get<double>();
        auto max_iter = json_solver["max_iter"].get<int>();
        auto directions = json_solver["directions"].get<std::array<int, 3>>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        LinearResponseSolver solver(dynamic, F_0, Phi, F_mat);
        solver.setHistory(kain);
        solver.setMaxIterations(max_iter);
        solver.setCanonical(not(localize));
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);

        F_0.setup(rsp_prec);
        for (int d = 0; d < 3; d++) {
            if (directions[d]) {
                F_1.perturbation() = h_1[d];
                X = orbital::param_copy(Phi);
                Y = orbital::param_copy(Phi);
                success = solver.optimize(omega, F_1, X, Y);
                if (success) {
                    const auto &json_prop = json_rsp["properties"].get<json>();
                    driver::calc_rsp_properties(json_prop, mol, d, omega);

                    auto final_orbs = json_rsp["final_orbitals"].get<std::string>();
                    auto write_orbs = json_rsp["write_orbitals"].get<bool>();
                    if (write_orbs) {
                        orbital::save_orbitals(X, final_orbs);
                        if (dynamic) orbital::save_orbitals(Y, final_orbs);
                    }
                }
                X.clear();
                Y.clear();
            }
            if (not(success)) break;
        }
        F_0.clear();
    }

    return success;
}

/** @brief Compute ground-state properties
 *
 * This function expects the "properties" subsection of the "scf_calculation"
 * input section, and will compute all properties which are present in this input.
 * This includes the diamagnetic contributions to the magnetic response properties.
 */
void driver::calc_scf_properties(const json &json_prop, Molecule &mol) {
    auto &nuclei = mol.getNuclei();
    auto &Phi = mol.getOrbitals();

    auto json_dipole = json_prop.find("dipole_moment");
    if (json_dipole != json_prop.end()) {
        Printer::printHeader(0, "Calculating dipole moment");
        auto prec = (*json_dipole)["setup_prec"].get<double>();
        auto r_O = (*json_dipole)["origin"].get<Coord<3>>();

        DipoleMoment &mu = mol.getDipoleMoment();

        Timer timer;
        H_E_dip h(r_O);
        h.setup(prec);
        mu.getOrigin() = r_O;
        mu.getNuclear() = h.trace(nuclei).real();
        mu.getElectronic() = h.trace(Phi).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }

    auto json_quadrupole = json_prop.find("quadrupole_moment");
    if (json_quadrupole != json_prop.end()) MSG_ERROR("Quadrupole moment not implemented");

    auto json_hyperpolarizability = json_prop.find("hyperpolarizability");
    if (json_hyperpolarizability != json_prop.end()) MSG_ERROR("Hyperpolarizability not implemented");

    auto json_gradient = json_prop.find("nuclear_gradient");
    if (json_gradient != json_prop.end()) {
        MSG_ERROR("Nuclear gradient not implemented");
        /*
        Printer::printHeader(0, "Calculating geometry derivatives");
        auto prec = json_input["grad_prec"].get<double>();
        auto smooth = json_input["grad_smooth"].get<double>();

        Timer timer;
        DoubleMatrix &nuc = mol.getGeometryDerivatives().getNuclear();
        DoubleMatrix &el = mol.getGeometryDerivatives().getElectronic();
        DoubleVector vecsum = DoubleVector::Zero(3);
        DoubleVector torque = DoubleVector::Zero(3);

        for (int k = 0; k < nuclei.size(); k++) {
            const Nucleus &nuc_k = nuclei[k];
            double Z_k = nuc_k.getCharge();
            const Coord<3> &R_k = nuc_k.getCoord();
            // Nuclear part
            for (int l = 0; l < nuclei.size(); l++) {
                if (l == k) continue;
                const Nucleus &nuc_l = nuclei[l];
                double Z_l = nuc_l.getCharge();
                const mrcpp::Coord<3> &R_l = nuc_l.getCoord();
                double R_kl = std::pow(math_utils::calc_distance(R_k, R_l), 3.0);
                nuc(k, 0) -= Z_k * Z_l * (R_k[0] - R_l[0]) / R_kl;
                nuc(k, 1) -= Z_k * Z_l * (R_k[1] - R_l[1]) / R_kl;
                nuc(k, 2) -= Z_k * Z_l * (R_k[2] - R_l[2]) / R_kl;
            }
            // Electronic part
            NuclearGradientOperator r_rm3(nuc_k, smooth);
            r_rm3.setup(prec);
            el.row(k) = r_rm3.trace(Phi).real();
            r_rm3.clear();

            // Total force
            vecsum += el.row(k);
            vecsum += nuc.row(k);

            // Total torque
            torque[0] += R_k[1] * (el(k, 2) + nuc(k, 2)) - R_k[2] * (el(k, 1) + nuc(k, 1));
            torque[1] += R_k[2] * (el(k, 0) + nuc(k, 0)) - R_k[0] * (el(k, 2) + nuc(k, 2));
            torque[2] += R_k[0] * (el(k, 1) + nuc(k, 1)) - R_k[1] * (el(k, 0) + nuc(k, 0));
        }
        println(0, "nuclear part    ");
        println(0, nuc);
        println(0, "electronic part ");
        println(0, el);
        println(0, "Total force acting on nuclei");
        println(0, vecsum.transpose());
        println(0, "Torque acting on nuclei");
        println(0, torque.transpose());
        timer.stop();
        Printer::printFooter(0, timer, 2);
    */
    }

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        Printer::printHeader(0, "Calculating diamagnetic magnetizability");
        auto prec = (*json_mag)["setup_prec"].get<double>();
        auto r_O = (*json_mag)["origin"].get<Coord<3>>();

        Magnetizability &khi = mol.getMagnetizability();

        Timer timer;
        H_BB_dia h(r_O);
        h.setup(prec);
        khi.getOrigin() = r_O;
        khi.getDiamagnetic() = -h.trace(Phi).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }

    auto json_nmr = json_prop.find("nmr_shielding");
    if (json_nmr != json_prop.end()) {
        Printer::printHeader(0, "Calculating diamagnetic NMR shielding");
        auto prec = (*json_nmr)["setup_prec"].get<double>();
        auto r_O = (*json_nmr)["origin"].get<Coord<3>>();
        auto nucleus_k = (*json_nmr)["nucleus_k"].get<std::vector<int>>();

        Timer timer;
        for (int k = 0; k < nucleus_k.size(); k++) {
            NMRShielding &sigma_k = mol.getNMRShielding(nucleus_k[k]);
            const auto &r_K = sigma_k.getNucleus().getCoord();

            H_BM_dia h(r_O, r_K);
            h.setup(prec);
            sigma_k.getOrigin() = r_O;
            sigma_k.getDiamagnetic() = h.trace(Phi).real();
            h.clear();
        }
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
}

/** @brief Compute linear response properties
 *
 * This function expects the "properties" subsection of the "rsp_calculations"
 * input section, and will compute all properties which are present in this input.
 */
void driver::calc_rsp_properties(const json &json_prop, Molecule &mol, int dir, double omega) {
    auto &Phi = mol.getOrbitals();
    auto &X = mol.getOrbitalsX();
    auto &Y = mol.getOrbitalsY();

    auto json_pol = json_prop.find("polarizability");
    if (json_pol != json_prop.end()) {
        Printer::printHeader(0, "Calculating polarizability");
        auto prec = (*json_pol)["setup_prec"].get<double>();
        auto r_O = (*json_pol)["origin"].get<Coord<3>>();

        Polarizability &alpha = mol.getPolarizability(omega);

        Timer timer;
        H_E_dip h(r_O);
        h.setup(prec);
        alpha.getTensor().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }

    auto json_mag = json_prop.find("magnetizability");
    if (json_mag != json_prop.end()) {
        Printer::printHeader(0, "Calculating paramagnetic magnetizability");
        auto prec = (*json_mag)["setup_prec"].get<double>();
        auto r_O = (*json_mag)["origin"].get<Coord<3>>();
        auto pert_diff = (*json_mag)["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);

        Magnetizability &khi = mol.getMagnetizability();

        Timer timer;
        H_B_dip h(D, r_O);
        h.setup(prec);
        khi.getParamagnetic().row(dir) = -h.trace(Phi, X, Y).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
}

/** @brief Build Fock operator based on input parameters
 *
 * This function expects the "fock_operator" subsection of input, and will
 * construct all operator which are present in this input. Option to set
 * perturbation order of the operators.
 */
void driver::build_fock_operator(const json &json_fock, Molecule &mol, FockOperator &F, int order) {
    auto &nuclei = mol.getNuclei();
    auto Phi_p = mol.getOrbitals_p();
    auto X_p = mol.getOrbitalsX_p();
    auto Y_p = mol.getOrbitalsY_p();

    ///////////////////////////////////////////////////////////
    //////////////////   Kinetic Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_kinetic = json_fock.find("kinetic_operator");
    if (json_kinetic != json_fock.end()) {
        auto kin_diff = (*json_kinetic)["derivative"].get<std::string>();
        auto D_p = driver::get_derivative(kin_diff);
        auto T_p = std::make_shared<KineticOperator>(D_p);
        F.getKineticOperator() = T_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Nuclear Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_nuclear = json_fock.find("nuclear_operator");
    if (json_nuclear != json_fock.end()) {
        auto proj_prec = (*json_nuclear)["proj_prec"].get<double>();
        auto smooth_prec = (*json_nuclear)["smooth_prec"].get<double>();
        auto shared_memory = (*json_nuclear)["shared_memory"].get<bool>();
        auto V_p = std::make_shared<NuclearOperator>(nuclei, proj_prec, smooth_prec, shared_memory);
        F.getNuclearOperator() = V_p;
    }
    ///////////////////////////////////////////////////////////
    //////////////////   Coulomb Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_coulomb = json_fock.find("coulomb_operator");
    if (json_coulomb != json_fock.end()) {
        auto poisson_prec = (*json_coulomb)["poisson_prec"].get<double>();
        auto shared_memory = (*json_coulomb)["shared_memory"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else if (order == 1) {
            auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, X_p, Y_p, shared_memory);
            F.getCoulombOperator() = J_p;
        } else {
            MSG_FATAL("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   Exchange Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    double exx = 1.0;
    auto json_exchange = json_fock.find("exchange_operator");
    if (json_exchange != json_fock.end()) {
        auto poisson_prec = (*json_exchange)["poisson_prec"].get<double>();
        auto screen_prec = (*json_exchange)["screen"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        if (order == 0) {
            auto K_p = std::make_shared<ExchangeOperator>(P_p, Phi_p, screen_prec);
            F.getExchangeOperator() = K_p;
        } else {
            MSG_FATAL("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    ////////////////////   XC Operator   //////////////////////
    ///////////////////////////////////////////////////////////
    auto json_xc = json_fock.find("xc_operator");
    if (json_xc != json_fock.end()) {
        auto grid_prec = (*json_xc)["grid_prec"].get<double>();
        auto shared_memory = (*json_xc)["shared_memory"].get<bool>();
        auto json_xcfunc = (*json_xc)["xc_functional"].get<json>();
        auto xc_spin = json_xcfunc["spin"].get<bool>();
        auto xc_gamma = json_xcfunc["gamma"].get<bool>();
        auto xc_cutoff = json_xcfunc["cutoff"].get<double>();
        auto xc_diff = json_xcfunc["derivative"].get<std::string>();
        auto xc_funcs = json_xcfunc["functionals"].get<json>();

        int xc_order = order + 1;
        auto xcfun_p = std::make_shared<XCFunctional>(*MRA, xc_spin);
        for (const auto &f : xc_funcs) {
            auto name = f["name"].get<std::string>();
            auto coef = f["coef"].get<double>();
            xcfun_p->setFunctional(name, coef);
        }
        xcfun_p->setUseGamma(xc_gamma);
        xcfun_p->setDensityCutoff(xc_cutoff);
        xcfun_p->evalSetup(xc_order);
        exx = xcfun_p->amountEXX();

        if (order == 0) {
            auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else if (order == 1) {
            auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, X_p, Y_p, shared_memory);
            F.getXCOperator() = XC_p;
        } else {
            MSG_FATAL("Invalid perturbation order");
        }
    }
    ///////////////////////////////////////////////////////////
    /////////////////   External Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_external = json_fock.find("external_operator");
    if (json_external != json_fock.end()) {
        auto field = (*json_external)["electric_field"].get<std::array<double, 3>>();
        auto r_O = (*json_external)["origin"].get<Coord<3>>();
        auto V_ext = std::make_shared<ElectricFieldOperator>(field, r_O);
        F.getExtOperator() = V_ext;
    }
    F.build(exx);
}

/** @brief Construct perturbation operator based on input keyword */
RankOneTensorOperator<3> driver::get_perturbation(const json &json_pert) {
    RankOneTensorOperator<3> h_1;
    auto pert_oper = json_pert["operator"].get<std::string>();
    if (pert_oper == "h_e_dip") {
        auto r_O = json_pert["origin"].get<mrcpp::Coord<3>>();
        h_1 = H_E_dip(r_O);
    }
    if (pert_oper == "h_b_dip") {
        auto r_O = json_pert["origin"].get<mrcpp::Coord<3>>();
        auto pert_diff = json_pert["derivative"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);
        h_1 = H_B_dip(D, r_O);
    }
    return h_1;
}

/** @brief Construct derivative operator based on input keyword */
DerivativeOperator_p driver::get_derivative(const std::string &name) {
    DerivativeOperator_p D = nullptr;
    if (name == "abgv_00") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.0, 0.0);
    } else if (name == "abgv_55") {
        D = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
        //    } else if (name == "ph") {
        //        D = std::make_shared<mrcpp::PHOperator<3>>(*MRA, 1);
        //    } else if (name == "bspline") {
        //        D = std::make_shared<mrcpp::BSOperator<3>>(*MRA, 1);
    } else {
        MSG_ERROR("Invalid derivative operator");
    }
    return D;
}

void driver::print_properties(const Molecule &mol) {
    println(0, "                                                            ");
    println(0, "************************************************************");
    println(0, "***                                                      ***");
    println(0, "***                Printing Properties                   ***");
    println(0, "***                                                      ***");
    println(0, "************************************************************");
    println(0, "                                                            ");
    println(0, "                                                            ");

    mol.printGeometry();
    mol.printProperties();
}

} // namespace mrchem

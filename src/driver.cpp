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
#include "scf_solver/OrbitalOptimizer.h"
#include "scf_solver/StaticResponseSolver.h"

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
void build_fock_operator(const json &input, Molecule &mol, FockOperator &F);
void build_perturbed_fock_operator(const json &input, Molecule &mol, FockOperator &F);

void calc_scf_properties(const json &input, Molecule &mol);
void calc_rsp_properties(const json &input, Molecule &mol);

RankOneTensorOperator<3> get_perturbation(const json &input);
DerivativeOperator_p get_derivative(const std::string &name);
} // namespace driver

void driver::init_molecule(const json &json_mol, Molecule &mol) {
    Printer::printHeader(0, "Molecule input");
    println(0, json_mol.dump(2));
    Printer::printSeparator(0, '=', 2);

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

bool driver::run_scf(const json &json_scf, Molecule &mol) {
    Printer::printHeader(0, "SCF input");
    println(0, json_scf.dump(2));
    Printer::printSeparator(0, '=', 2);

    auto success = true;
    auto scf_prec = json_scf["scf_prec"].get<double>();

    const auto &json_fock = json_scf["fock_operator"].get<json>();
    FockOperator F;
    driver::build_fock_operator(json_fock, mol, F);

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
        auto canonical = json_solver["canonical"].get<bool>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        OrbitalOptimizer solver;
        solver.setHistory(kain);
        solver.setMaxIterations(max_iter);
        solver.setRotation(rotation);
        solver.setCanonical(canonical);
        solver.setOrbitalPrec(start_prec, final_prec);
        solver.setThreshold(orbital_thrs, property_thrs);
        success = solver.optimize(F, Phi, F_mat);
    } else {
        F.setup(scf_prec);
        F_mat = F(Phi, Phi);
        F.clear();
    }

    // Run EnergyOptimizer if present in input JSON
    auto energy_solver_it = json_scf.find("energy_solver");
    if (energy_solver_it != json_scf.end()) {
        const auto &json_solver = *energy_solver_it;
        auto max_iter = json_solver["max_iter"].get<int>();
        auto canonical = json_solver["canonical"].get<bool>();
        auto start_prec = json_solver["start_prec"].get<double>();
        auto final_prec = json_solver["final_prec"].get<double>();
        auto orbital_thrs = json_solver["orbital_thrs"].get<double>();
        auto property_thrs = json_solver["property_thrs"].get<double>();

        EnergyOptimizer solver;
        solver.setMaxIterations(max_iter);
        solver.setRotation(1);
        solver.setCanonical(canonical);
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

bool driver::run_rsp(const json &json_rsp, Molecule &mol) {
    NOT_IMPLEMENTED_ABORT;
}

void driver::calc_scf_properties(const json &json_prop, Molecule &mol) {
    Nuclei &nuclei = mol.getNuclei();
    OrbitalVector &Phi = mol.getOrbitals();

    auto json_dipole = json_prop.find("dipole_moment");
    if (json_dipole != json_prop.end()) {
        Printer::printHeader(0, "Calculating dipole moment");
        auto prec = (*json_dipole)["setup_prec"].get<double>();
        auto r_O = (*json_dipole)["origin"].get<Coord<3>>();

        Timer timer;
        DipoleMoment &mu = mol.getDipoleMoment();

        H_E_dip h(r_O);
        h.setup(prec);
        mu.getNuclear() = h.trace(nuclei).real();
        mu.getElectronic() = h.trace(Phi).real();
        h.clear();
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }

    auto json_quadrupole = json_prop.find("quadrupole_moment");
    if (json_quadrupole != json_prop.end()) MSG_ERROR("Quadrupole moment not implemented");

    auto json_polarizability = json_prop.find("polarizability");
    if (json_polarizability != json_prop.end()) MSG_ERROR("Polarizability not implemented");

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

        Timer timer;
        Magnetizability &eta = mol.getMagnetizability();

        H_BB_dia h(r_O);
        h.setup(prec);
        eta.getDiamagnetic() = -h.trace(Phi).real();
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
            sigma_k.getDiamagnetic() = h.trace(Phi).real();
            h.clear();
        }
        timer.stop();
        Printer::printFooter(0, timer, 2);
    }
}

void calc_rsp_properties(const json &json_input, Molecule &mol) {
    NOT_IMPLEMENTED_ABORT;
}

void driver::build_fock_operator(const json &json_fock, Molecule &mol, FockOperator &F) {
    auto &nuclei = mol.getNuclei();
    auto Phi_p = mol.getOrbitals_p();

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
        auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, shared_memory);
        F.getCoulombOperator() = J_p;
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
        auto K_p = std::make_shared<ExchangeOperator>(P_p, Phi_p, screen_prec);
        F.getExchangeOperator() = K_p;
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

        int xc_order = 1;
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

        auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, shared_memory);
        F.getXCOperator() = XC_p;
    }
    F.build(exx);
}

void driver::build_perturbed_fock_operator(const json &json_fock, Molecule &mol, FockOperator &F) {
    auto Phi_p = mol.getOrbitals_p();
    auto X_p = mol.getOrbitalsX_p();
    auto Y_p = mol.getOrbitalsY_p();

    ///////////////////////////////////////////////////////////
    //////////////////   Coulomb Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    auto json_coulomb = json_fock.find("coulomb_operator");
    if (json_coulomb != json_fock.end()) {
        auto poisson_prec = (*json_coulomb)["poisson_prec"].get<double>();
        auto shared_memory = (*json_coulomb)["shared_memory"].get<bool>();
        auto P_p = std::make_shared<PoissonOperator>(*MRA, poisson_prec);
        auto J_p = std::make_shared<CoulombOperator>(P_p, Phi_p, X_p, Y_p, shared_memory);
        F.getCoulombOperator() = J_p;
    }
    ///////////////////////////////////////////////////////////
    /////////////////   Exchange Operator   ///////////////////
    ///////////////////////////////////////////////////////////
    double exx = 1.0;
    auto json_exchange = json_fock.find("exchange_operator");
    if (json_exchange != json_fock.end()) NOT_IMPLEMENTED_ABORT;
    ///////////////////////////////////////////////////////////
    ////////////////////   XC Operator   //////////////////////
    ///////////////////////////////////////////////////////////
    auto json_xc = json_fock.find("xc_operator");
    if (json_xc != json_fock.end()) {
        auto grid_prec = (*json_xc)["grid_prec"].get<double>();
        auto shared_memory = (*json_xc)["shared_memory"].get<bool>();
        auto json_xcfunc = (*json_xc)["xcfunctional"].get<json>();
        auto xc_spin = json_xcfunc["spin"].get<bool>();
        auto xc_gamma = json_xcfunc["gamma"].get<bool>();
        auto xc_cutoff = json_xcfunc["cutoff"].get<double>();
        auto xc_diff = json_xcfunc["derivative"].get<std::string>();
        auto xc_funcs = json_xcfunc["functionals"].get<json>();

        int xc_order = 2;
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

        auto XC_p = std::make_shared<XCOperator>(xcfun_p, Phi_p, X_p, Y_p, shared_memory);
        F.getXCOperator() = XC_p;
    }
    F.build(exx);
}

RankOneTensorOperator<3> driver::get_perturbation(const json &json_pert) {
    RankOneTensorOperator<3> h_1;
    auto pert_oper = json_pert["pert_oper"].get<std::string>();
    if (pert_oper == "h_e_dip") {
        auto r_O = json_pert["gauge_origin"].get<mrcpp::Coord<3>>();
        h_1 = H_E_dip(r_O);
    }
    if (pert_oper == "h_b_dip") {
        auto r_O = json_pert["gauge_origin"].get<mrcpp::Coord<3>>();
        auto pert_diff = json_pert["pert_diff"].get<std::string>();
        auto D = driver::get_derivative(pert_diff);
        h_1 = H_B_dip(D, r_O);
    }
    return h_1;
}

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

} // namespace mrchem

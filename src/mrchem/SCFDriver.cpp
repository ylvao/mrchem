#include "SCFDriver.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include "eigen_disable_warnings.h"

#include "TelePrompter.h"

#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"

#include "MWAdder.h"
#include "MWMultiplier.h"
#include "MWProjector.h"
#include "GridGenerator.h"
#include "DerivativeOperator.h"
#include "PoissonOperator.h"

#include "CoreHamiltonian.h"
//#include "Hartree.h"
//#include "HartreeFock.h"
//#include "DFT.h"

//#include "GroundStateSolver.h"
//#include "LinearResponseSolver.h"
//#include "HelmholtzOperatorSet.h"
//#include "KAIN.h"
//#include "DIIS.h"

#include "Molecule.h"
#include "OrbitalVector.h"
//#include "Orbital.h"
//#include "Density.h"
#include "Potential.h"

#include "InitialGuessProjector.h"

//#include "SCFEnergy.h"
//#include "DipoleMoment.h"
//#include "QuadrupoleMoment.h"
//#include "Polarizability.h"
//#include "OpticalRotation.h"
//#include "Magnetizability.h"
//#include "NMRShielding.h"
//#include "SpinSpinCoupling.h"

#include "KineticOperator.h"
#include "NuclearFunction.h"
//#include "CoulombPotential.h"
//#include "CoulombHessian.h"
//#include "ExchangePotential.h"
//#include "ExchangeHessian.h"
//#include "XCFunctional.h"
//#include "XCPotential.h"
//#include "XCHessian.h"
//#include "DipoleOperator.h"
//#include "QuadrupoleOperator.h"
//#include "AngularMomentumOperator.h"
//#include "MomentumOperator.h"

#include "mrchem.h"

using namespace std;
using namespace Eigen;

SCFDriver::SCFDriver(Getkw &input) {
    order = input.get<int>("order");
    rel_prec = input.get<double>("rel_prec");
    max_depth = input.get<int>("max_depth");

    scale = input.get<int>("World.scale");
    center_of_mass = input.get<bool>("World.center_of_mass");
    boxes = input.getIntVec("World.boxes");
    corner = input.getIntVec("World.corner");
    gauge = input.getDblVec("World.gauge_origin");

    run_ground_state = input.get<bool>("Properties.ground_state");
    run_dipole_moment = input.get<bool>("Properties.dipole_moment");
    run_quadrupole_moment = input.get<bool>("Properties.quadrupole_moment");
    run_polarizability = input.get<bool>("Properties.polarizability");
    run_optrot_electric = input.get<bool>("Properties.optrot_electric");
    run_optrot_magnetic = input.get<bool>("Properties.optrot_magnetic");
    run_magnetizability = input.get<bool>("Properties.magnetizability");
    run_nmr_shielding = input.get<bool>("Properties.nmr_shielding");
    run_spin_spin = input.get<bool>("Properties.spin_spin_coupling");
    nmr_nuclei = input.getIntVec("Properties.nmr_nuclei");
    spin_spin_k = input.getIntVec("Properties.spin_spin_k");
    spin_spin_l = input.getIntVec("Properties.spin_spin_l");
    frequencies = input.getDblVec("Properties.frequencies");
    velocity_gauge = input.get<bool>("Properties.velocity_gauge");

    mol_charge = input.get<int>("Molecule.charge");
    mol_multiplicity = input.get<int>("Molecule.multiplicity");
    mol_coords = input.getData("Molecule.coords");

    wf_restricted = input.get<bool>("WaveFunction.restricted");
    wf_method = input.get<string>("WaveFunction.method");

    if (wf_method == "DFT") {
        dft_spin = input.get<bool>("DFT.spin");
        dft_x_fac = input.get<double>("DFT.exact_exchange");
        dft_cutoff = input.getDblVec("DFT.density_cutoff");
        dft_func_coefs = input.getDblVec("DFT.func_coefs");
        dft_func_names = input.getData("DFT.functionals");
    }

    scf_start = input.get<string>("SCF.initial_guess");
    scf_history = input.get<int>("SCF.history");
    scf_max_iter = input.get<int>("SCF.max_iter");
    scf_rotation = input.get<int>("SCF.rotation");
    scf_localize = input.get<bool>("SCF.localize");
    scf_write_orbitals = input.get<bool>("SCF.write_orbitals");
    scf_orbital_thrs = input.get<double>("SCF.orbital_thrs");
    scf_property_thrs = input.get<double>("SCF.property_thrs");
    scf_lambda_thrs = input.get<double>("SCF.lambda_thrs");
    scf_orbital_prec = input.getDblVec("SCF.orbital_prec");

    rsp_dir = input.getIntVec("Response.directions");
    rsp_start = input.get<string>("Response.initial_guess");
    rsp_history = input.get<int>("Response.history");
    rsp_max_iter = input.get<int>("Response.max_iter");
    rsp_localize = input.get<bool>("Response.localize");
    rsp_write_orbitals = input.get<bool>("Response.write_orbitals");
    rsp_orbital_thrs = input.get<double>("Response.orbital_thrs");
    rsp_property_thrs = input.get<double>("Response.property_thrs");
    rsp_orbital_prec = input.getDblVec("Response.orbital_prec");

    file_start_orbitals = input.get<string>("Files.start_orbitals");
    file_final_orbitals = input.get<string>("Files.final_orbitals");
    file_start_x_orbs = input.get<string>("Files.start_x_orbs");
    file_final_x_orbs = input.get<string>("Files.final_x_orbs");
    file_start_y_orbs = input.get<string>("Files.start_y_orbs");
    file_final_y_orbs = input.get<string>("Files.final_y_orbs");
    file_basis_set = input.get<string>("Files.basis_set");
    file_dens_mat = input.get<string>("Files.dens_mat");
    file_fock_mat = input.get<string>("Files.fock_mat");
    file_energy_vec = input.get<string>("Files.energy_vec");
    file_mo_mat_a = input.get<string>("Files.mo_mat_a");
    file_mo_mat_b = input.get<string>("Files.mo_mat_b");

    run_el_field_rsp = false;
    run_mag_field_rsp = false;
    run_mag_moment_rsp = false;
    if (run_nmr_shielding) run_mag_field_rsp = true;
    if (run_magnetizability) run_mag_field_rsp = true;
    if (run_optrot_magnetic) run_mag_field_rsp = true;
    if (run_polarizability) run_el_field_rsp = true;
    if (run_optrot_electric) run_el_field_rsp = true;
    if (run_spin_spin) run_mag_moment_rsp = true;

    r_O[0] = 0.0;
    r_O[1] = 0.0;
    r_O[2] = 0.0;

    helmholtz = 0;
    scf_kain = 0;
    rsp_kain_x = 0;
    rsp_kain_y = 0;

    molecule = 0;
    phi = 0;
    x_phi = 0;
    y_phi = 0;

    P = 0;
    T = 0;
    V = 0;
    J = 0;
    K = 0;
    XC = 0;

    dJ = 0;
    dK = 0;
    dXC = 0;

    f_mat = 0;
    f_oper = 0;
    df_oper = 0;

    xcfun_1 = 0;
    xcfun_2 = 0;
}

bool SCFDriver::sanityCheck() const {
    if (wf_restricted and wf_method == "DFT" and dft_spin) {
        MSG_ERROR("Restricted spin DFT not implemented");
        return false;
    }
    if (wf_restricted and mol_multiplicity != 1) {
        MSG_ERROR("Restricted open-shell not implemented");
        return false;
    }
    if (not wf_restricted and wf_method == "HF") {
        MSG_ERROR("Unrestricted HF not implemented");
        return false;
    }
    if (not wf_restricted and (run_mag_field_rsp or run_el_field_rsp)) {
        MSG_ERROR("Unrestricted response not implemented");
        return false;
    }
    return true;
}

void SCFDriver::setup() {
    // Setting up MRA
    NodeIndex<3> idx(scale, corner.data());
    BoundingBox<3> world(idx, boxes.data());
    InterpolatingBasis basis(order);
    MRA = new MultiResolutionAnalysis<3>(world, basis, max_depth);

    // Setting up MW operators
    add = new MWAdder<3>(*MRA);
    mult = new MWMultiplier<3>(*MRA);
    Q = new MWProjector<3>(*MRA, rel_prec);
    G = new GridGenerator<3>(*MRA);
    D = new DerivativeOperator<3>(*MRA, rel_prec, rel_prec/10.0);
    P = new PoissonOperator(*MRA, rel_prec, rel_prec/10.0);

    // Setting up molecule
    molecule = new Molecule(mol_coords, mol_charge);
    int nEl = molecule->getNElectrons();
    nuclei = &molecule->getNuclei();
    phi = new OrbitalVector(nEl, mol_multiplicity, wf_restricted);

    // Defining gauge origin
    const double *COM = molecule->getCenterOfMass();
    if (center_of_mass) {
        r_O[0] = COM[0];
        r_O[1] = COM[1];
        r_O[2] = COM[2];
    } else {
        r_O[0] = gauge[0];
        r_O[1] = gauge[1];
        r_O[2] = gauge[2];
    }

    // Setup properties
    if (nmr_nuclei[0] < 0) {
        nmr_nuclei.clear();
        for (int i = 0; i < molecule->getNNuclei(); i++) {
            nmr_nuclei.push_back(i);
        }
    }
    if (spin_spin_k[0] < 0) {
        spin_spin_k.clear();
        for (int k = 0; k < molecule->getNNuclei(); k++) {
            spin_spin_k.push_back(k);
        }
    }
    if (spin_spin_l[0] < 0) {
        spin_spin_l.clear();
        for (int l = 0; l < molecule->getNNuclei(); l++) {
            spin_spin_l.push_back(l);
        }
    }
    if (run_dipole_moment) {
        molecule->initDipoleMoment(r_O);
    }
    if (run_quadrupole_moment) {
        molecule->initQuadrupoleMoment(r_O);
    }
    for (int i = 0; i < frequencies.size(); i++) {
        if (run_polarizability) {
            molecule->initPolarizability(frequencies[i], r_O, velocity_gauge);
        }
        if (run_optrot_electric or run_optrot_magnetic) {
            molecule->initOpticalRotation(frequencies[i], r_O, velocity_gauge);
        }
    }
    if (run_magnetizability) {
        molecule->initMagnetizability(0.0, r_O);
    }
    if (run_nmr_shielding) {
        for (int k = 0; k < nmr_nuclei.size(); k++) {
            int K = nmr_nuclei[k];
            molecule->initNMRShielding(K, r_O);
        }
    }
    if (run_spin_spin) {
        for (int k = 0; k < spin_spin_k.size(); k++) {
            int K = spin_spin_k[k];
            for (int l = 0; l < spin_spin_l.size(); l++) {
                int L = spin_spin_l[l];
                molecule->initSpinSpinCoupling(K, L);
            }
        }
    }

    // Setting up Fock matrix
    int nOrbs = phi->size();
    f_mat = new MatrixXd;
    *f_mat = MatrixXd::Zero(nOrbs, nOrbs);

    // Setting up SCF
//    helmholtz = new HelmholtzOperatorSet(scf_lambda_thrs);
//    if (scf_history > 0) scf_kain = new KAIN(scf_history);
//    if (rsp_history > 0) rsp_kain_x = new KAIN(rsp_history);
//    if (rsp_history > 0) rsp_kain_y = new KAIN(rsp_history);

    Timer nuc_timer;
    nuc_timer.restart();
    TelePrompter::printHeader(0, "Setting up nuclear potential");
    NuclearFunction nuc_func(*nuclei, rel_prec);
    nuc_tree = (*Q)(nuc_func);
    TelePrompter::printFooter(0, nuc_timer, 2);

    // Setting up Fock operator
    T = new KineticOperator(*D);
    V = new Potential(*add, *mult, nuc_tree, 0);
//    rho = new DensityOperator(*phi);
//    J = new CoulombPotential(*P, *phi);

    if (wf_method == "Core") {
        f_oper = new CoreHamiltonian(*T, *V);
//    } else if (wf_method == "Hartree") {
//        f_oper = new Hartree(*T, *V, *J);
//    } else if (wf_method == "HF") {
//        K = new ExchangePotential(*P, *phi);
//        f_oper = new HartreeFock(*T, *V, *J, *K);
//    } else if (wf_method == "DFT") {
//        xcfun_1 = new XCFunctional(dft_spin, 1);
//        xcfun_1->setDensityCutoff(dft_cutoff[0]);
//        for (int i = 0; i < dft_func_names.size(); i++) {
//            xcfun_1->setFunctional(dft_func_names[i], dft_func_coefs[i]);
//        }
//        XC = new XCPotential(*xcfun_1, *phi);
//        if (dft_x_fac > MachineZero) {
//            K = new ExchangePotential(*P, *phi, dft_x_fac);
//        }
//        f_oper = new DFT(*T, *V, *J, *XC, K);
    } else {
        MSG_ERROR("Invalid method");
    }
}

void SCFDriver::clear() {
    if (MRA != 0) delete MRA;

    if (molecule != 0) delete molecule;
    if (phi != 0) delete phi;
    if (nuc_tree != 0) delete nuc_tree;

//    if (helmholtz != 0) delete helmholtz;
//    if (scf_kain != 0) delete scf_kain;
//    if (rsp_kain_x != 0) delete rsp_kain_x;
//    if (rsp_kain_y != 0) delete rsp_kain_y;

    if (add != 0) delete add;
    if (mult != 0) delete mult;
    if (Q != 0) delete Q;
    if (G != 0) delete G;
    if (D != 0) delete D;
    if (P != 0) delete P;
//    if (T != 0) delete T;
//    if (J != 0) delete J;
//    if (V != 0) delete V;
//    if (K != 0) delete K;
//    if (XC != 0) delete XC;
    if (f_mat != 0) delete f_mat;
    if (f_oper != 0) delete f_oper;
//    if (xcfun_1 != 0) delete xcfun_1;
}

GroundStateSolver* SCFDriver::setupInitialGuessSolver() {
    NOT_IMPLEMENTED_ABORT;
//    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");

//    GroundStateSolver *gss = new GroundStateSolver(*helmholtz);
//    gss->setMaxIterations(-1);
//    gss->setRotation(-1);
//    gss->setThreshold(1.0e-1, 1.0);
//    gss->setOrbitalPrec(1.0e-3, -1.0);
//    return gss;
}

GroundStateSolver* SCFDriver::setupGroundStateSolver() {
    NOT_IMPLEMENTED_ABORT;
//    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");

//    GroundStateSolver *gss = new GroundStateSolver(*helmholtz, scf_kain);
//    gss->setMaxIterations(scf_max_iter);
//    gss->setRotation(scf_rotation);
//    gss->setThreshold(scf_orbital_thrs, scf_property_thrs);
//    gss->setOrbitalPrec(scf_orbital_prec[0], scf_orbital_prec[1]);
//    return gss;
}

LinearResponseSolver* SCFDriver::setupLinearResponseSolver(bool dynamic) {
    NOT_IMPLEMENTED_ABORT;
//    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");
//    LinearResponseSolver *lrs = 0;
//    if (dynamic) {
//        lrs = new LinearResponseSolver(*helmholtz, rsp_kain_x, rsp_kain_y);
//    } else {
//        lrs = new LinearResponseSolver(*helmholtz, rsp_kain_x);
//    }
//    lrs->setMaxIterations(rsp_max_iter);
//    lrs->setThreshold(rsp_orbital_thrs, rsp_property_thrs);
//    lrs->setOrbitalPrec(rsp_orbital_prec[0], rsp_orbital_prec[1]);
//    lrs->setupUnperturbedSystem(*f_oper, *f_mat, *orbitals);
//    return lrs;
}

void SCFDriver::setupPerturbedOrbitals(bool dynamic) {
    NOT_IMPLEMENTED_ABORT;
//    if (phi == 0) MSG_ERROR("Orbitals not initialized");

//    x_phi = new OrbitalVector("xOrbs", *phi);
//    x_phi->initialize(*phi);
//    if (dynamic) {
//        y_phi = new OrbitalVector("yOrbs", *phi);
//        y_phi->initialize(*phi);
//    } else {
//        y_phi = x_phi;
//    }
}

void SCFDriver::clearPerturbedOrbitals(bool dynamic) {
    NOT_IMPLEMENTED_ABORT;
//    if (x_phi != 0) delete x_phi;
//    if (dynamic) {
//        if (y_phi != 0) delete y_phi;
//    }
//    x_phi = 0;
//    y_phi = 0;
}

void SCFDriver::setupPerturbedOperators(QMOperator &dH, bool dynamic) {
    NOT_IMPLEMENTED_ABORT;
//    if (P == 0) MSG_ERROR("Poission operator not initialized");
//    if (phi == 0) MSG_ERROR("Orbitals not initialized");
//    if (x_phi == 0) MSG_ERROR("X orbitals not initialized");
//    if (y_phi == 0) MSG_ERROR("Y orbitals not initialized");

//    double xFac = 0.0;
//    if (wf_method == "HF") {
//        xFac = 1.0;
//    } else if (wf_method == "DFT") {
//        xFac = dft_x_fac;
//    }
//    if (xFac > MachineZero) {
//        dK = new ExchangeHessian(*P, *phi, x_phi, y_phi, xFac);
//    }

//    bool imaginary = dH.isImaginary();
//    if (not imaginary or dynamic) {
//        dJ = new CoulombHessian(*P, *phi, x_phi, y_phi);
//        if (wf_method == "DFT") {
//            xcfun_2 = new XCFunctional(dft_spin, 2);
//            xcfun_2->setDensityCutoff(dft_cutoff[1]);
//            for (int i = 0; i < dft_func_names.size(); i++) {
//                xcfun_2->setFunctional(dft_func_names[i], dft_func_coefs[i]);
//            }
//            if (xcfun_2 == 0) MSG_ERROR("xcfun not initialized");
//            dXC = new XCHessian(*xcfun_2, *phi, x_phi, y_phi);
//        }
//    }

//    df_oper = new FockOperator(0, 0, dJ, dK, dXC);
//    df_oper->addPerturbationOperator(dH);
}

void SCFDriver::clearPerturbedOperators() {
    NOT_IMPLEMENTED_ABORT;
//    if (dJ != 0) delete dJ;
//    if (dK != 0) delete dK;
//    if (dXC != 0) delete dXC;
//    if (df_oper != 0) delete df_oper;
//    if (xcfun_2 != 0) delete xcfun_2;

//    dJ = 0;
//    dK = 0;
//    dXC = 0;
//    df_oper = 0;
//    xcfun_2 = 0;
}

void SCFDriver::run() {
    bool converged = true;
    if (not sanityCheck()) {
        return;
    }
    setupInitialGroundState();
    if (run_ground_state) {
        converged = runGroundState();
    } else {
        f_oper->setup();
        *f_mat = (*f_oper)(*phi, *phi);
        f_oper->clear();
    }
//    calcGroundStateProperties();
//    if (converged) {
//        if (run_el_field_rsp or run_mag_field_rsp) {
//            if (rsp_localize) {
//                phi->localize(-1.0, f_mat);
//            } else {
//                phi->diagonalize(-1.0, f_mat);
//            }
//        }
//        if (run_el_field_rsp) {
//            for (int i = 0; i < frequencies.size(); i++) {
//                runElectricFieldResponse(frequencies[i]);
//            }
//        }
//        if (run_mag_field_rsp) {
//            for (int i = 0; i < frequencies.size(); i++) {
//                runMagneticFieldResponse(frequencies[i]);
//            }
//        }
//        if (run_mag_moment_rsp) {
//            for (int i = 0; i < spin_spin_l.size(); i++) {
//                int L = spin_spin_l[i];
//                runMagneticMomentResponse("PSO", L);
//                runMagneticMomentResponse("SD", L);
//                runMagneticMomentResponse("FC", L);
//            }
//        }
//    }
    TelePrompter::printHeader(0, "Fock matrix");
    println(0, *f_mat);
    TelePrompter::printSeparator(0, '=', 2);
    molecule->printGeometry();
//    molecule->printProperties();
}

bool SCFDriver::runInitialGuess(FockOperator &oper, MatrixXd &F, OrbitalVector &orbs) {
    NOT_IMPLEMENTED_ABORT;
//    if (phi == 0) MSG_ERROR("Orbitals not initialized");

//    GroundStateSolver *gss = setupInitialGuessSolver();
//    gss->setup(oper, F, orbs);
//    bool converged = gss->optimize();
//    gss->clear();

//    delete gss;
//    return converged;
}

bool SCFDriver::runGroundState() {
    NOT_IMPLEMENTED_ABORT;
//    if (phi == 0) MSG_ERROR("Orbitals not initialized");
//    if (f_oper == 0) MSG_ERROR("Fock operator not initialized");
//    if (f_mat == 0) MSG_ERROR("Fock matrix not initialized");

//    GroundStateSolver *gss = setupGroundStateSolver();
//    gss->setup(*f_oper, *f_mat, *phi);
//    bool converged = gss->optimize();
//    gss->clear();

//    if (scf_write_orbitals) {
//        orbitals->writeOrbitals(file_final_orbitals);
//    }

//    delete gss;
//    return converged;
}

void SCFDriver::runElectricFieldResponse(double omega) {
    NOT_IMPLEMENTED_ABORT;
//    bool dynamic = false;
//    if (omega > MachineZero) dynamic = true;

//    LinearResponseSolver *lrs = setupLinearResponseSolver(dynamic);

//    bool converged = false;
//    for (int d = 0; d < 3; d++) {
//        if (rsp_dir[d] < 1) continue;
//        QMOperator *h = 0;
//        if (velocity_gauge) {
//            h = new MomentumOperator(d);
//        } else {
//            h = new DipoleOperator(d, r_O);
//        }

//        setupPerturbedOrbitals(dynamic);
//        setupInitialResponse(*h, d, dynamic, omega);
//        setupPerturbedOperators(*h, dynamic);

//        if (dynamic) {
//            lrs->setup(omega, *df_oper, *x_phi, *y_phi);
//        } else {
//            lrs->setup(*df_oper, *x_phi);
//        }
//        converged = lrs->optimize();
//        lrs->clear();

//        if (converged) calcElectricFieldProperties(d, omega);

//        if (rsp_write_orbitals) {
//            stringstream x_file;
//            x_file << file_final_x_orbs << "_" << d;
//            x_phi->writeOrbitals(x_file.str());
//            if (dynamic) {
//                stringstream y_file;
//                y_file << file_final_y_orbs << "_" << d;
//                y_phi->writeOrbitals(y_file.str());
//            }
//        }

//        clearPerturbedOperators();
//        clearPerturbedOrbitals(dynamic);
//        delete h;

//        if (not converged) break;
//    }

//    if (lrs != 0) lrs->clearUnperturbedSystem();
//    if (lrs != 0) delete lrs;
}

void SCFDriver::runMagneticFieldResponse(double omega) {
    NOT_IMPLEMENTED_ABORT;
//    bool dynamic = false;
//    if (omega > MachineZero) dynamic = true;

//    LinearResponseSolver *lrs = setupLinearResponseSolver(dynamic);

//    bool converged = false;
//    for (int d = 0; d < 3; d++) {
//        if (rsp_dir[d] < 1) continue;
//        AngularMomentumOperator h(d, r_O);

//        setupPerturbedOrbitals(dynamic);
//        setupInitialResponse(h, d, dynamic, omega);
//        setupPerturbedOperators(h, dynamic);

//        if (dynamic) {
//            lrs->setup(omega, *df_oper, *x_phi, *y_phi);
//        } else {
//            lrs->setup(*df_oper, *x_phi);
//        }
//        converged = lrs->optimize();
//        lrs->clear();

//        if (converged) calcMagneticFieldProperties(d, omega);

//        if (rsp_write_orbitals) {
//            stringstream x_file;
//            x_file << file_final_x_orbs << "_" << d;
//            x_phi->writeOrbitals(x_file.str());
//            if (dynamic) {
//                stringstream y_file;
//                y_file << file_final_y_orbs << "_" << d;
//                y_phi->writeOrbitals(y_file.str());
//            }
//        }

//        clearPerturbedOperators();
//        clearPerturbedOrbitals(dynamic);

//        if (not converged) break;
//    }

//    if (lrs != 0) lrs->clearUnperturbedSystem();
//    if (lrs != 0) delete lrs;
}

void SCFDriver::runMagneticMomentResponse(const string &type, int L) {
    NOT_IMPLEMENTED_ABORT;
//    double omega = 0.0;
//    bool dynamic = false;
//    LinearResponseSolver *lrs = setupLinearResponseSolver(dynamic);
//    const double *r_L = (*nuclei)[L]->getCoord();

//    bool converged = false;
//    for (int d = 0; d < 3; d++) {
//        if (rsp_dir[d] < 1) continue;
//        QMOperator *h = 0;
//        if (type == "PSO") h = new PSOOperator(d, r_L);
//        if (type == "SD") NOT_IMPLEMENTED_ABORT;
//        if (type == "FC") NOT_IMPLEMENTED_ABORT;
//        if (h == 0) MSG_ERROR("Invalid type");

//        setupPerturbedOrbitals(dynamic);
//        setupInitialResponse(*h, d, dynamic, omega);
//        setupPerturbedOperators(*h, dynamic);

//        lrs->setup(*df_oper, *x_orbs);
//        converged = lrs->optimize();
//        lrs->clear();

//        if (converged) calcMagneticMomentProperties(type, d, omega);

//        if (rsp_write_orbitals) {
//            stringstream x_file;
//            x_file << file_final_x_orbs << "_" << d;
//            x_phi->writeOrbitals(x_file.str());
//            if (dynamic) {
//                stringstream y_file;
//                y_file << file_final_y_orbs << "_" << d;
//                y_phi->writeOrbitals(y_file.str());
//            }
//        }

//        clearPerturbedOperators();
//        clearPerturbedOrbitals(dynamic);
//        if (h != 0) delete h;

//        if (not converged) break;
//    }

//    if (lrs != 0) lrs->clearUnperturbedSystem();
//    if (lrs != 0) delete lrs;
}

void SCFDriver::calcGroundStateProperties() {
    NOT_IMPLEMENTED_ABORT;
//    f_oper->setup();
//    SCFEnergy &energy = molecule->getSCFEnergy();
//    energy.compute(*nuclei);
//    energy.compute(*f_oper, *f_mat, *phi);
//    f_oper->clear();

//    if (run_dipole_moment) {
//        DipoleMoment &dipole = molecule->getDipoleMoment();
//        dipole.compute(*nuclei);
//        dipole.compute(*phi);
//    }
//    if (run_quadrupole_moment) {
//        QuadrupoleMoment &quadrupole = molecule->getQuadrupoleMoment();
//        quadrupole.compute(*nuclei);
//        quadrupole.compute(*phi);
//    }
//    if (run_magnetizability) {
//        Magnetizability &magnetizability = molecule->getMagnetizability();
//        magnetizability.computeDiamagnetic(*phi);
//    }
//    if (run_nmr_shielding) {
//        for (int k = 0; k < nmr_nuclei.size(); k++) {
//            int K = nmr_nuclei[k];
//            NMRShielding &shielding = molecule->getNMRShielding(K);
//            shielding.computeDiamagnetic(*phi);
//        }
//    }
//    if (run_spin_spin) {
//        for (int k = 0; k < spin_spin_k.size(); k++) {
//            int K = spin_spin_k[k];
//            for (int l = 0; l < spin_spin_l.size(); l++) {
//                int L = spin_spin_l[l];
//                SpinSpinCoupling &ssc = molecule->getSpinSpinCoupling(K, L);
//                ssc.computeDiamagnetic(*phi);
//            }
//        }
//    }
}

void SCFDriver::calcElectricFieldProperties(int d, double omega) {
    NOT_IMPLEMENTED_ABORT;
//    if (run_polarizability) {
//        Polarizability &polarizability = molecule->getPolarizability(omega);
//        polarizability.compute(d, *phi, *x_phi, *y_phi);
//    }
//    if (run_optrot_electric) {
//        OpticalRotation &opt_rot = molecule->getOpticalRotation(omega);
//        opt_rot.compute("E", d, *phi, *x_phi, *y_phi);
//    }
}

void SCFDriver::calcMagneticFieldProperties(int d, double omega) {
    NOT_IMPLEMENTED_ABORT;
//    if (run_magnetizability) {
//        Magnetizability &magnetizability = molecule->getMagnetizability();
//        magnetizability.computeParamagnetic("B", d, *phi, *x_phi, *y_phi);
//    }
//    if (run_nmr_shielding) {
//        for (int k = 0; k < nmr_nuclei.size(); k++) {
//            int K = nmr_nuclei[k];
//            NMRShielding &shielding = molecule->getNMRShielding(K);
//            shielding.computeParamagnetic("B", d, *phi, *x_phi, *y_phi);
//        }
//    }
//    if (run_optrot_magnetic) {
//        OpticalRotation &opt_rot = molecule->getOpticalRotation(omega);
//        opt_rot.compute("B", d, *phi, *x_phi, *y_phi);
//    }
}

void SCFDriver::calcMagneticMomentProperties(const string &type, int d, int L) {
    NOT_IMPLEMENTED_ABORT;
//    if (run_spin_spin) {
//        for (int k = 0; k < nmr_nuclei.size(); k++) {
//            int K = nmr_nuclei[k];
//            SpinSpinCoupling &ssc = molecule->getSpinSpinCoupling(K, L);
//            ssc.computeParamagnetic(type, d, *phi, *x_phi, *y_phi);
//        }
//    }
}

void SCFDriver::setupInitialGroundState() {
    // Reading initial guess
    if (scf_start == "none") {
        NOT_IMPLEMENTED_ABORT;
//        OrbitalVector initOrbs("initOrbs");
//        initOrbs.initialize(*nuclei);
//        MatrixXd S = initOrbs.calcOverlapMatrix();
//        println(0, endl << S << endl);
//        CoreHamiltonian h(*T, *V);
//        MatrixXd F = h(initOrbs, initOrbs);
//        println(0, endl << F << endl);
//        initOrbs.diagonalize(-1.0, &F);
//        println(0, endl << F << endl);
//        phi->readOrbitals(initOrbs);
//        runInitialGuess(*f_oper, F, *phi);
    } else if (scf_start == "gto") {
        InitialGuessProjector IGP(*MRA, rel_prec);
        if (wf_restricted) {
            IGP(*phi, file_basis_set, file_mo_mat_a);
        } else {
            IGP(*phi, file_basis_set, file_mo_mat_a, file_mo_mat_b);
        }
    } else if (scf_start == "mw") {
        NOT_IMPLEMENTED_ABORT;
//        phi->readOrbitals(file_start_orbitals);
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
}

void SCFDriver::setupInitialResponse(QMOperator &h, int d,
                                     bool dynamic, double omega) {
    NOT_IMPLEMENTED_ABORT;
    // Reading initial guess
//    if (rsp_start == "none") {
//    } else if (rsp_start == "gto") {
//        int nOcc = phi->size();
//        OrbitalVector *mol_orbs = new OrbitalVector("virtuals");
//        mol_orbs->readVirtuals(file_basis_set, file_mo_mat_a, nOcc);
//        int nVirt = mol_orbs->size();

//        if (rsp_localize) phi->diagonalize(-1.0, f_mat);
//        VectorXd epsilon = MathUtils::readVectorFile(file_energy_vec);

//        for (int i = 0; i < nOcc; i++) {
//            vector<double> x_ai;
//            vector<double> y_ai;
//            vector<FunctionTree<3> *> t_a;
//            Orbital &phi_i = phi->getOrbital(i);
//            for (int a = 0; a < nVirt; a++) {
//                Orbital &phi_a = mol_orbs->getOrbital(a);
//                double h_ai = -h(phi_a, phi_i);
//                if (fabs(h_ai) > MachineZero) {
//                    x_ai.push_back(h_ai/(epsilon(nOcc + a) - (*f_mat)(i,i) - omega));
//                    y_ai.push_back(h_ai/(epsilon(nOcc + a) - (*f_mat)(i,i) + omega));
//                    t_a.push_back(&phi_a);
//                }
//            }
//            if (t_a.size() > 0) {
//                Orbital &x_i = x_phi->getOrbital(i);
//                x_i.add(x_ai, t_a, 0);
//                printout(0, "xOrb    " << setw(3) << i);
//                println(0, " squareNorm: " << setw(36) << x_i.getSquareNorm());
//                if (dynamic) {
//                    Orbital &y_i = y_phi->getOrbital(i);
//                    y_i.add(y_ai, t_a, 0);
//                }
//            }
//        }
//        if (rsp_localize) {
//            MatrixXd U = phi->localize(-1.0, f_mat);
//            x_phi->rotate(-1.0, U);
//            if (dynamic) y_phi->rotate(-1.0, U);
//        }
//        delete mol_orbs;
//    } else if (rsp_start == "mw") {
//        stringstream x_file, y_file;
//        x_file << file_final_x_orbs << "_" << d;
//        x_phi->readOrbitals(x_file.str());
//        if (dynamic) {
//            y_file << file_final_y_orbs << "_" << d;
//            y_phi->readOrbitals(y_file.str());
//        }
//    } else {
//        NOT_IMPLEMENTED_ABORT;
//    }
}

/*
void SCFDriver::runMagneticallyInducedCurrents() {
    int oldPrec = TelePrompter::setPrecision(5);
    printout(0, "\n------------------------------------------------------------\n");
    printout(0, "\nGauge Origin:   ");
    printout(0, setw(14)<<rsp_gauge[0]<<setw(14)<<rsp_gauge[1]<<setw(14)<<rsp_gauge[2]<<endl);
    printout(0, "\n------------------------------------------------------------\n");
    TelePrompter::setPrecision(oldPrec);

    initResponse(frequencies);

    DerivativeOperator del(0.0, 0.0);
    OrbitalVector &orbs = molecule.getOrbitals();

// calculate perturbed orbitals
    OrbitalVector *pOrbs[3];
    pOrbs[0] = new OrbitalVector("xOrbs", orbs);
    pOrbs[1] = new OrbitalVector("yOrbs", orbs);
    pOrbs[2] = new OrbitalVector("zOrbs", orbs);
    calcMagneticResponseX(0, r_O, orbs, *pOrbs[0]);
    calcMagneticResponseX(1, r_O, orbs, *pOrbs[1]);
    calcMagneticResponseX(2, r_O, orbs, *pOrbs[2]);

// setup current tensors
    FunctionTree<3> **j_tot[3];
    FunctionTree<3> **j_dia[3];
    FunctionTree<3> **j_para1[3];
    FunctionTree<3> **j_para2[3];
    for (int i = 0; i < 3; i++) {
        j_tot[i] = new FunctionTree<3>*[3];
        j_dia[i] = new FunctionTree<3>*[3];
        j_para1[i] = new FunctionTree<3>*[3];
        j_para2[i] = new FunctionTree<3>*[3];
        for (int j = 0; j < 3; j++) {
            j_tot[i][j] = new FunctionTree<3>;
            j_dia[i][j] = new FunctionTree<3>;
            j_para1[i][j] = new FunctionTree<3>;
            j_para2[i][j] = new FunctionTree<3>;
        }
    }

// setup x-, y-, z-functions for the diamagnetic contribution
    FunctionTree<3> *rTree[3];
    PositionFunction<3> rFunc(r_O);
    for (int d = 0; d < 3; d++) {
        rTree[d] = new FunctionTree<3>;
        rFunc.setDir(d);
        rTree[d]->projectFunction(rFunc);
    }

    Density *rho = new Density(-1.0, orbs, orbs);
// compute the diamagnetic contribution to the current
    j_dia[0][1]->mult(-0.5, *rTree[2], 1.0, *rho);
    j_dia[0][2]->mult( 0.5, *rTree[1], 1.0, *rho);

    j_dia[1][0]->mult( 0.5, *rTree[2], 1.0, *rho);
    j_dia[1][2]->mult(-0.5, *rTree[0], 1.0, *rho);

    j_dia[2][0]->mult(-0.5, *rTree[1], 1.0, *rho);
    j_dia[2][1]->mult( 0.5, *rTree[0], 1.0, *rho);

    delete rho;
    delete rTree[0];
    delete rTree[1];
    delete rTree[2];

// compute the paramagnetic contribution to the current
    for (int n = 0; n < orbs.size(); n++) {
        Orbital &orb = orbs.getOrbital(n);
        for (int j = 0; j < 3; j++) {
            Orbital dOrb(orb);
            del.setApplyDir(j);
            del.apply(dOrb, orb);

            for (int i = 0; i < 3; i++) {
                Orbital &pOrb = pOrbs[i]->getOrbital(n);
                Orbital dpOrb(orb);
                del.apply(dpOrb, pOrb);

                FunctionTree<3> orbPart1;
                FunctionTree<3> orbPart2;
                orbPart1.mult( 2.0, orb, 1.0, dpOrb);
                orbPart2.mult(-2.0, pOrb, 1.0, dOrb);

                j_para1[i][j]->addInPlace(1.0, orbPart1);
                j_para2[i][j]->addInPlace(1.0, orbPart2);
            }
        }
    }

    // Evaluating the contributions to the current
    const double *origin = BoundingBox<3>::getWorldBox().getOrigin();
    double r[3] = {0.1 - origin[0],
                   0.2 - origin[1],
                   0.3 - origin[2] };

    Matrix3d int_dia;
    Matrix3d int_para1;
    Matrix3d int_para2;
    Matrix3d eval_dia;
    Matrix3d eval_para1;
    Matrix3d eval_para2;
    for(int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            int_dia(i,j) = j_dia[i][j]->integrate();
            int_para1(i,j) = j_para1[i][j]->integrate();
            int_para2(i,j) = j_para2[i][j]->integrate();
            eval_dia(i,j) = j_dia[i][j]->evalf(r);
            eval_para1(i,j) = j_para1[i][j]->evalf(r);
            eval_para2(i,j) = j_para2[i][j]->evalf(r);
        }
    }

    println(0,"");
    println(0, "integral diamagnetic:");
    println(0, int_dia);
    println(0,"");
    println(0, "integral para1:");
    println(0, int_para1);
    println(0,"");
    println(0, "integral para2:");
    println(0, int_para2);
    println(0,"");
    println(0, "integral paramagnetic:");
    println(0, int_para1 + int_para2);
    println(0,"");
    println(0, "integral total");
    println(0, int_dia + int_para1 + int_para2);

    println(0,"");
    println(0,"");
    println(0, "diamagnetic:");
    println(0, eval_dia);

    println(0,"");
    println(0, "para1:");
    println(0, eval_para1);

    println(0,"");
    println(0, "para2:");
    println(0, eval_para2);

    println(0,"");
    println(0, "paramagnetic:");
    println(0, eval_para1 + eval_para2);

    println(0,"");
    println(0, "total:");
    println(0, eval_dia + eval_para1 + eval_para2);

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            delete j_tot[i][j];
            delete j_dia[i][j];
            delete j_para1[i][j];
            delete j_para2[i][j];
        }
        delete[] j_tot[i];
        delete[] j_dia[i];
        delete[] j_para1[i];
        delete[] j_para2[i];
    }

    delete pOrbs[0];
    delete pOrbs[1];
    delete pOrbs[2];
}
*/


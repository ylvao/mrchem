#include "SCFDriver.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Eigenvalues>

#include "mrchem.h"
#include "TelePrompter.h"
#include "MathUtils.h"
#include "eigen_disable_warnings.h"

#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"

#include "CoreHamiltonian.h"
#include "Hartree.h"
#include "HartreeFock.h"
#include "DFT.h"

#include "OrbitalOptimizer.h"
#include "EnergyOptimizer.h"
#include "HelmholtzOperatorSet.h"
#include "KAIN.h"

#include "Molecule.h"
#include "OrbitalVector.h"

#include "OrbitalProjector.h"

#include "SCFEnergy.h"
#include "DipoleMoment.h"

#include "DipoleOperator.h"

#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombPotential.h"
#include "ExchangePotential.h"
#include "XCPotential.h"
#include "XCFunctional.h"

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

    mol_charge = input.get<int>("Molecule.charge");
    mol_multiplicity = input.get<int>("Molecule.multiplicity");
    mol_coords = input.getData("Molecule.coords");

    wf_restricted = input.get<bool>("WaveFunction.restricted");
    wf_method = input.get<string>("WaveFunction.method");

    if (wf_method == "DFT") {
        dft_spin = input.get<bool>("DFT.spin");
        dft_x_fac = input.get<double>("DFT.exact_exchange");
        dft_cutoff = input.get<double>("DFT.density_cutoff");
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

    file_start_orbitals = input.get<string>("Files.start_orbitals");
    file_final_orbitals = input.get<string>("Files.final_orbitals");
    file_basis_set = input.get<string>("Files.basis_set");
    file_dens_mat = input.get<string>("Files.dens_mat");
    file_fock_mat = input.get<string>("Files.fock_mat");
    file_energy_vec = input.get<string>("Files.energy_vec");
    file_mo_mat_a = input.get<string>("Files.mo_mat_a");
    file_mo_mat_b = input.get<string>("Files.mo_mat_b");

    r_O[0] = 0.0;
    r_O[1] = 0.0;
    r_O[2] = 0.0;

    helmholtz = 0;
    scf_kain = 0;

    molecule = 0;
    phi = 0;

    T = 0;
    V = 0;
    J = 0;
    K = 0;
    XC = 0;
    fock = 0;

    J_np1 = 0;
    K_np1 = 0;
    XC_np1 = 0;
    fock_np1 = 0;

    xcfun = 0;
}

bool SCFDriver::sanityCheck() const {
    if (not wf_restricted) {
        MSG_ERROR("Unrestricted SCF not implemented");
        return false;
    }
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
    return true;
}

void SCFDriver::setup() {
    // Setting up MRA
    NodeIndex<3> idx(scale, corner.data());
    BoundingBox<3> world(idx, boxes.data());
    InterpolatingBasis basis(order);
    MRA = new MultiResolutionAnalysis<3>(world, basis, max_depth);

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

    if (run_dipole_moment) {
        molecule->initDipoleMoment(r_O);
    }

    // Setting up SCF
    helmholtz = new HelmholtzOperatorSet(rel_prec, *MRA, scf_lambda_thrs);
    if (scf_history > 0) scf_kain = new KAIN(*MRA, scf_history);

    // Setting up Fock operator
    T = new KineticOperator(rel_prec, *MRA);
    V = new NuclearPotential(rel_prec, *MRA, *nuclei);

    if (wf_method == "Core") {
        fock = new CoreHamiltonian(*MRA, *T, *V);
    } else if (wf_method == "Hartree") {
        J = new CoulombPotential(rel_prec, *MRA, *phi);
        fock = new Hartree(*MRA, *T, *V, *J);
    } else if (wf_method == "HF") {
        J = new CoulombPotential(rel_prec, *MRA, *phi);
        K = new ExchangePotential(rel_prec, *MRA, *phi);
        fock = new HartreeFock(*MRA, *T, *V, *J, *K);
    } else if (wf_method == "DFT") {
        J = new CoulombPotential(rel_prec, *MRA, *phi);
        xcfun = new XCFunctional(dft_spin, dft_cutoff);
        for (int i = 0; i < dft_func_names.size(); i++) {
            xcfun->setFunctional(dft_func_names[i], dft_func_coefs[i]);
        }
        XC = new XCPotential(rel_prec, *MRA, *xcfun, *phi);
        if (dft_x_fac > MachineZero) {
            K = new ExchangePotential(rel_prec, *MRA, *phi, dft_x_fac);
        }
        fock = new DFT(*MRA, *T, *V, *J, *XC, 0);
    } else {
        MSG_ERROR("Invalid method");
    }
}

void SCFDriver::clear() {
    if (MRA != 0) delete MRA;

    if (molecule != 0) delete molecule;
    if (phi != 0) delete phi;

    if (helmholtz != 0) delete helmholtz;
    if (scf_kain != 0) delete scf_kain;

    if (T != 0) delete T;
    if (V != 0) delete V;
    if (J != 0) delete J;
    if (K != 0) delete K;
    if (XC != 0) delete XC;
    if (fock != 0) delete fock;
    if (xcfun != 0) delete xcfun;
}

/** Setup n+1 Fock operator for energy optimization */
void SCFDriver::setup_np1() {
    phi_np1 = new OrbitalVector(*phi);

    if (wf_method == "Core") {
    } else if (wf_method == "Hartree") {
        J_np1 = new CoulombPotential(rel_prec, *MRA, *phi_np1);
    } else if (wf_method == "HF") {
        J_np1 = new CoulombPotential(rel_prec, *MRA, *phi_np1);
        K_np1 = new ExchangePotential(rel_prec, *MRA, *phi_np1);
    } else if (wf_method == "DFT") {
        J_np1 = new CoulombPotential(rel_prec, *MRA, *phi_np1);
        XC_np1 = new XCPotential(rel_prec, *MRA, *xcfun, *phi_np1);
        if (dft_x_fac > MachineZero) {
            K_np1 = new ExchangePotential(rel_prec, *MRA, *phi_np1, dft_x_fac);
        }
    } else {
        MSG_ERROR("Invalid method");
    }

    fock_np1 = new FockOperator(*MRA, 0, V, J_np1, K_np1, XC_np1);
}

void SCFDriver::clear_np1() {
    if (phi_np1 != 0) delete phi_np1;
    if (J_np1 != 0) delete J_np1;
    if (K_np1 != 0) delete K_np1;
    if (XC_np1 != 0) delete XC_np1;
    if (fock_np1 != 0) delete fock_np1;
}

void SCFDriver::setupInitialGroundState() {
    // Reading initial guess
    if (scf_start == "none") {
        // Project minimal basis set of hydrogen orbitals
        OrbitalProjector OP(*MRA, rel_prec);
        OrbitalVector *tmp = OP(*nuclei);

        // Compute orthonormalization matrix
        MatrixXd S = tmp->calcOverlapMatrix().real();
        println(0, endl << S << endl);
        MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S, -1.0/2.0);

        // Compute core Hamiltonian matrix
        CoreHamiltonian h(*MRA, *T, *V);
        h.setup(rel_prec);
        MatrixXd f_mat = h(*tmp, *tmp);
        h.clear();
        println(0, endl << f_mat << endl);

        // Diagonalize core Hamiltonian matrix
        MatrixXd M = MathUtils::diagonalizeHermitianMatrix(f_mat);
        MatrixXd U = M.transpose()*S_m12;
        println(0, endl << f_mat << endl);

        // Rotate n lowest energy orbitals of U*tmp into phi
        OrbitalAdder add(*MRA, rel_prec);
        add.rotate(*phi, U, *tmp);
        delete tmp;
    } else if (scf_start == "gto") {
        OrbitalProjector OP(*MRA, rel_prec);
        if (wf_restricted) {
            OP(*phi, file_basis_set, file_mo_mat_a);
        } else {
            OP(*phi, file_basis_set, file_mo_mat_a, file_mo_mat_b);
        }
    } else if (scf_start == "mw") {
        NOT_IMPLEMENTED_ABORT;
    } else {
        NOT_IMPLEMENTED_ABORT;
    }
}

OrbitalOptimizer* SCFDriver::setupOrbitalOptimizer() {
    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");

    OrbitalOptimizer *optimizer = new OrbitalOptimizer(*MRA, *helmholtz, scf_kain);
    optimizer->setMaxIterations(scf_max_iter);
    optimizer->setRotation(scf_rotation);
    optimizer->setThreshold(scf_orbital_thrs, scf_property_thrs);
    optimizer->setOrbitalPrec(scf_orbital_prec[0], scf_orbital_prec[1]);

    return optimizer;
}

EnergyOptimizer* SCFDriver::setupEnergyOptimizer() {
    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");

    EnergyOptimizer *optimizer = new EnergyOptimizer(*MRA, *helmholtz);
    optimizer->setMaxIterations(scf_max_iter);
    optimizer->setRotation(0);
    optimizer->setThreshold(scf_orbital_thrs, scf_property_thrs);
    optimizer->setOrbitalPrec(scf_orbital_prec[0], scf_orbital_prec[1]);

    return optimizer;
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
        fock->setup(rel_prec);
        F = (*fock)(*phi, *phi);
        fock->clear();
    }
    calcGroundStateProperties();

    printEigenvalues(F, *phi);
    molecule->printGeometry();
    molecule->printProperties();
}

bool SCFDriver::runGroundState() {
    if (phi == 0) MSG_ERROR("Orbitals not initialized");
    if (fock == 0) MSG_ERROR("Fock operator not initialized");

    bool converged = false;
    {   // Optimize orbitals
        OrbitalOptimizer *solver = setupOrbitalOptimizer();
        solver->setup(*fock, *phi, F);
        converged = solver->optimize();
        solver->clear();
        delete solver;
    }

        // Optimize energy
    if (scf_property_thrs > 0.0) {
        setup_np1();

        EnergyOptimizer *solver = setupEnergyOptimizer();
        solver->setup(*fock, *phi, F, *fock_np1, *phi_np1);
        converged = solver->optimize();
        solver->clear();

        clear_np1();
        delete solver;
    }

    if (scf_write_orbitals) {
        NOT_IMPLEMENTED_ABORT;
//        phi->writeOrbitals(file_final_orbitals);
    }

    return converged;
}

void SCFDriver::calcGroundStateProperties() {
    SCFEnergy &energy = molecule->getSCFEnergy();
    fock->setup(rel_prec);
    energy.compute(*nuclei);
    energy.compute(*fock, F, *phi);
    fock->clear();

    if (run_dipole_moment) {
        DipoleMoment &dipole = molecule->getDipoleMoment();
        for (int d = 0; d < 3; d++) {
            DipoleOperator mu_d(*MRA, d, r_O[d]);
            mu_d.setup(rel_prec);
            dipole.compute(d, mu_d, *nuclei);
            dipole.compute(d, mu_d, *phi);
            mu_d.clear();
        }
    }
}

void SCFDriver::printEigenvalues(MatrixXd &f_mat, OrbitalVector &orbs) {
    int oldprec = TelePrompter::setPrecision(5);
    TelePrompter::printHeader(0, "Fock matrix");
    println(0, f_mat);
    TelePrompter::printSeparator(0, '=', 2);

    TelePrompter::printHeader(0, "Orbital energies");
    println(0, "    n  spin  occ                            epsilon  ");
    TelePrompter::printSeparator(0, '-');
    SelfAdjointEigenSolver<MatrixXd> es(f_mat.cols());
    es.compute(f_mat);

    TelePrompter::setPrecision(15);
    VectorXd epsilon = es.eigenvalues();
    for (int i = 0; i < epsilon.size(); i++) {
        Orbital &orb = orbs.getOrbital(i);
        printout(0, setw(5) << i);
        printout(0, setw(5) << orb.printSpin());
        printout(0, setw(5) << orb.getOccupancy());
        printout(0, setw(44) << epsilon(i) << endl);
    }
    TelePrompter::printSeparator(0, '=', 2);
    TelePrompter::setPrecision(oldprec);
}

#include <iostream>
#include <iomanip>
#include <fstream>
#include <Eigen/Eigenvalues>

#include "eigen_disable_warnings.h"
#include "mrchem.h"
#include "SCFDriver.h"
#include "TelePrompter.h"
#include "MathUtils.h"

#include "InterpolatingBasis.h"
#include "MultiResolutionAnalysis.h"

#include "CoreHamiltonian.h"
#include "Hartree.h"
#include "HartreeFock.h"
#include "DFT.h"

#include "HelmholtzOperatorSet.h"
#include "OrbitalOptimizer.h"
#include "EnergyOptimizer.h"
#include "KAIN.h"

#include "Molecule.h"
#include "OrbitalVector.h"
#include "OrbitalProjector.h"

#include "SCFEnergy.h"
#include "DipoleMoment.h"
#include "Magnetizability.h"
#include "NMRShielding.h"
#include "SpinSpinCoupling.h"

#include "DipoleOperator.h"
#include "DMOperator.h"
#include "DSOperator.h"
#include "DSOOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombPotential.h"
#include "ExchangePotential.h"
#include "XCPotential.h"
#include "XCFunctional.h"

using namespace std;
using namespace Eigen;

SCFDriver::SCFDriver(Getkw &input) {
    rel_prec = input.get<double>("rel_prec");

    gauge = input.getDblVec("World.gauge_origin");
    center_of_mass = input.get<bool>("World.center_of_mass");

    calc_total_energy = input.get<bool>("Properties.total_energy");
    calc_dipole_moment = input.get<bool>("Properties.dipole_moment");
    calc_quadrupole_moment = input.get<bool>("Properties.quadrupole_moment");
    calc_magnetizability = input.get<bool>("Properties.magnetizability");
    calc_nmr_shielding = input.get<bool>("Properties.nmr_shielding");
    calc_hyperfine_coupling = input.get<bool>("Properties.hyperfine_coupling");
    calc_spin_spin_coupling = input.get<bool>("Properties.spin_spin_coupling");
    calc_polarizability = input.get<bool>("Properties.polarizability");
    calc_hyperpolarizability = input.get<bool>("Properties.hyperpolarizability");
    calc_optical_rotation = input.get<bool>("Properties.optical_rotation");

    nmr_perturbation = input.get<string>("NMRShielding.perturbation");
    nmr_nucleus_k = input.getIntVec("NMRShielding.nucleus_k");
    hfcc_nucleus_k = input.getIntVec("HyperfineCoupling.nucleus_k");
    sscc_nucleus_k = input.getIntVec("SpinSpinCoupling.nucleus_k");
    sscc_nucleus_l = input.getIntVec("SpinSpinCoupling.nucleus_l");
    pol_velocity = input.get<bool>("Polarizability.velocity");
    pol_frequency = input.getDblVec("Polarizability.frequency");
    optrot_velocity = input.get<bool>("OpticalRotation.velocity");
    optrot_frequency = input.getDblVec("OpticalRotation.frequency");
    optrot_perturbation = input.get<string>("OpticalRotation.perturbation");

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

    scf_run = input.get<bool>("SCF.run");
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
    if (calc_quadrupole_moment) {
        MSG_ERROR("Quadrupole moment not implemented");
        return false;
    }
    if (calc_polarizability) {
        MSG_ERROR("Polarizability not implemented");
        return false;
    }
    if (calc_hyperpolarizability) {
        MSG_ERROR("Hyperpolarizability not implemented");
        return false;
    }
    if (calc_optical_rotation) {
        MSG_ERROR("Optical rotation not implemented");
        return false;
    }
    if (calc_magnetizability) {
        MSG_ERROR("Only diamagnetic magnetizability atm");
    }
    if (calc_nmr_shielding) {
        MSG_ERROR("Only diamagnetic NMR shielding atm");
    }
    if (calc_spin_spin_coupling) {
        MSG_ERROR("Only diamagnetic spin-spin coupling atm");
    }
    if (calc_hyperfine_coupling) {
        MSG_ERROR("Hyperfine coupling not implemented");
        return false;
    }

    return true;
}

void SCFDriver::setup() {
    // Setting up molecule
    molecule = new Molecule(mol_coords, mol_charge);
    int nEl = molecule->getNElectrons();
    nuclei = &molecule->getNuclei();

    // Setting up empty orbitals
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

    // Setting up properties
    if (nmr_nucleus_k[0] < 0) {
        nmr_nucleus_k.clear();
        for (int k = 0; k < molecule->getNNuclei(); k++) {
            nmr_nucleus_k.push_back(k);
        }
    }
    if (sscc_nucleus_k[0] < 0) {
        sscc_nucleus_k.clear();
        for (int k = 0; k < molecule->getNNuclei(); k++) {
            sscc_nucleus_k.push_back(k);
        }
    }
    if (sscc_nucleus_l[0] < 0) {
        sscc_nucleus_l.clear();
        for (int l = 0; l < molecule->getNNuclei(); l++) {
            sscc_nucleus_l.push_back(l);
        }
    }
    if (hfcc_nucleus_k[0] < 0) {
        hfcc_nucleus_k.clear();
        for (int k = 0; k < molecule->getNNuclei(); k++) {
            hfcc_nucleus_k.push_back(k);
        }
    }

    if (calc_dipole_moment) molecule->initDipoleMoment();
    if (calc_quadrupole_moment) molecule->initQuadrupoleMoment();
    if (calc_magnetizability) molecule->initMagnetizability();
    if (calc_nmr_shielding) {
        for (int k = 0; k < nmr_nucleus_k.size(); k++) {
            int K = nmr_nucleus_k[k];
            molecule->initNMRShielding(K);
        }
    }
    if (calc_hyperfine_coupling) {
        for (int k = 0; k < hfcc_nucleus_k.size(); k++) {
            int K = hfcc_nucleus_k[k];
            molecule->initHyperfineCoupling(K);
        }
    }
    if (calc_spin_spin_coupling) {
        for (int k = 0; k < sscc_nucleus_k.size(); k++) {
            int K = sscc_nucleus_k[k];
            for (int l = 0; l < sscc_nucleus_l.size(); l++) {
                int L = sscc_nucleus_l[l];
                if (K != L) molecule->initSpinSpinCoupling(K, L);
            }
        }
    }
    if (calc_polarizability) {
        for (int i = 0; i < pol_frequency.size(); i++) {
            molecule->initPolarizability(pol_frequency[i]);
        }
    }
    if (calc_optical_rotation) {
        for (int i = 0; i < optrot_frequency.size(); i++) {
            molecule->initOpticalRotation(optrot_frequency[i]);
        }
    }

    // Setting up SCF
    helmholtz = new HelmholtzOperatorSet(rel_prec, scf_lambda_thrs);
    if (scf_history > 0) scf_kain = new KAIN(scf_history);

    // Setting up Fock operator
    T = new KineticOperator(rel_prec);
    V = new NuclearPotential(rel_prec, *nuclei);

    if (wf_method == "Core") {
        fock = new CoreHamiltonian(*T, *V);
    } else if (wf_method == "Hartree") {
        J = new CoulombPotential(rel_prec, *phi);
        fock = new Hartree(*T, *V, *J);
    } else if (wf_method == "HF") {
        J = new CoulombPotential(rel_prec, *phi);
        K = new ExchangePotential(rel_prec, *phi);
        fock = new HartreeFock(*T, *V, *J, *K);
    } else if (wf_method == "DFT") {
        J = new CoulombPotential(rel_prec, *phi);
        xcfun = new XCFunctional(dft_spin, dft_cutoff);
        for (int i = 0; i < dft_func_names.size(); i++) {
            xcfun->setFunctional(dft_func_names[i], dft_func_coefs[i]);
        }
        XC = new XCPotential(rel_prec, *xcfun, *phi);
        if (dft_x_fac > MachineZero) {
            K = new ExchangePotential(rel_prec, *phi, dft_x_fac);
        }
        fock = new DFT(*T, *V, *J, *XC, 0);
    } else {
        MSG_ERROR("Invalid method");
    }
}

void SCFDriver::clear() {
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
        J_np1 = new CoulombPotential(rel_prec, *phi_np1);
    } else if (wf_method == "HF") {
        J_np1 = new CoulombPotential(rel_prec, *phi_np1);
        K_np1 = new ExchangePotential(rel_prec, *phi_np1);
    } else if (wf_method == "DFT") {
        J_np1 = new CoulombPotential(rel_prec, *phi_np1);
        XC_np1 = new XCPotential(rel_prec, *xcfun, *phi_np1);
        if (dft_x_fac > MachineZero) {
            K_np1 = new ExchangePotential(rel_prec, *phi_np1, dft_x_fac);
        }
    } else {
        MSG_ERROR("Invalid method");
    }

    fock_np1 = new FockOperator(0, V, J_np1, K_np1, XC_np1);
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
        OrbitalProjector OP(rel_prec);
        OrbitalVector *tmp = OP(*nuclei);

        // Compute orthonormalization matrix
        MatrixXd S = tmp->calcOverlapMatrix().real();
        MatrixXd S_m12 = MathUtils::hermitianMatrixPow(S, -1.0/2.0);

        // Compute core Hamiltonian matrix
        CoreHamiltonian h(*T, *V);
        h.setup(rel_prec);
        MatrixXd f_mat = h(*tmp, *tmp);
        h.clear();

        // Diagonalize core Hamiltonian matrix
        MatrixXd M = MathUtils::diagonalizeHermitianMatrix(f_mat);
        MatrixXd U = M.transpose()*S_m12;

        // Rotate n lowest energy orbitals of U*tmp into phi
        OrbitalAdder add(rel_prec);
        add.rotate(*phi, U, *tmp);
        delete tmp;
    } else if (scf_start == "gto") {
        OrbitalProjector OP(rel_prec);
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

    OrbitalOptimizer *optimizer = new OrbitalOptimizer(*helmholtz, scf_kain);
    optimizer->setMaxIterations(scf_max_iter);
    optimizer->setRotation(scf_rotation);
    optimizer->setThreshold(scf_orbital_thrs, scf_property_thrs);
    optimizer->setOrbitalPrec(scf_orbital_prec[0], scf_orbital_prec[1]);

    return optimizer;
}

EnergyOptimizer* SCFDriver::setupEnergyOptimizer() {
    if (helmholtz == 0) MSG_ERROR("Helmholtz operators not initialized");

    EnergyOptimizer *optimizer = new EnergyOptimizer(*helmholtz);
    optimizer->setMaxIterations(scf_max_iter);
    if (scf_localize) {
        optimizer->setRotation(1);
    } else {
        optimizer->setRotation(-1);
    }
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
    if (scf_run) {
        converged = runGroundState();
    } else {
        fock->setup(rel_prec);
        F = (*fock)(*phi, *phi);
        fock->clear();
    }
    calcGroundStateProperties();

    printEigenvalues(*phi, F);
    molecule->printGeometry();
    molecule->printProperties();
}

bool SCFDriver::runGroundState() {
    if (phi == 0) MSG_ERROR("Orbitals not initialized");
    if (fock == 0) MSG_ERROR("Fock operator not initialized");

    bool converged = false;

    // Optimize orbitals
    if (scf_orbital_thrs > 0.0) {
        OrbitalOptimizer *solver = setupOrbitalOptimizer();
        solver->setup(*fock, *phi, F);
        converged = solver->optimize();
        solver->clear();
        delete solver;
    } else {
        fock->setup(rel_prec);
        F = (*fock)(*phi, *phi);
        fock->clear();
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
    }

    return converged;
}

void SCFDriver::calcGroundStateProperties() {
    if (calc_total_energy) {
        SCFEnergy &energy = molecule->getSCFEnergy();
        fock->setup(rel_prec);
        energy.compute(*nuclei);
        energy.compute(*fock, *phi, F);
        fock->clear();
    }

    if (calc_dipole_moment) {
        Vector3d &nuc = molecule->getDipoleMoment().getNuclear();
        Vector3d &el = molecule->getDipoleMoment().getElectronic();
        for (int i = 0; i < 3; i++) {
            DipoleOperator mu_i(i, r_O);
            mu_i.setup(rel_prec);
            nuc(i) = mu_i.trace(*nuclei);
            el(i) = mu_i.trace(*phi);
            mu_i.clear();
        }
    }
    if (calc_magnetizability) {
        Matrix3d &dia = molecule->getMagnetizability().getDiamagnetic();
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                DMOperator h_BB(i, j, r_O);
                h_BB.setup(rel_prec);
                dia(i,j) = h_BB.trace(*phi);
                h_BB.clear();
            }
        }
    }
    if (calc_nmr_shielding) {
        for (int k = 0; k < nmr_nucleus_k.size(); k++) {
            int K = nmr_nucleus_k[k];
            NMRShielding &nmr = molecule->getNMRShielding(K);
            Matrix3d &dia = nmr.getDiamagnetic();
            const double *r_K = nmr.getNucleus().getCoord();
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    DSOperator h_BM(i, j, r_O, r_K);
                    h_BM.setup(rel_prec);
                    dia(i,j) = h_BM.trace(*phi);
                    h_BM.clear();
                }
            }
        }
    }
    if (calc_spin_spin_coupling) {
        for (int k = 0; k < sscc_nucleus_k.size(); k++) {
            int K = sscc_nucleus_k[k];
            for (int l = 0; l < sscc_nucleus_l.size(); l++) {
                int L = sscc_nucleus_l[l];
                if (K == L) continue;
                SpinSpinCoupling &sscc = molecule->getSpinSpinCoupling(K, L);
                Matrix3d &dia = sscc.getDiamagnetic();
                const double *r_K = sscc.getNucleusK().getCoord();
                const double *r_L = sscc.getNucleusL().getCoord();
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        DSOOperator h_MM(i, j, r_K, r_L);
                        h_MM.setup(rel_prec);
                        dia(i,j) = h_MM.trace(*phi);
                        h_MM.clear();
                    }
                }
            }
        }
    }
}

void SCFDriver::printEigenvalues(OrbitalVector &orbs, MatrixXd &f_mat) {
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

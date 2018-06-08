#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "initial_guess/sad.h"
#include "initial_guess/gto.h"
#include "initial_guess/core.h"

#include "Molecule.h"

#include "NuclearOperator.h"
#include "KineticOperator.h"
#include "CoulombOperator.h"
#include "XCOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace sad {
void project_atomic_densities(double prec, const Molecule &mol, DensityVector &rho_atomic);
} //namespace sad
} //namespace initial_guess

OrbitalVector initial_guess::sad::setup(double prec,
                                            const Molecule &mol,
                                            bool restricted,
                                            int zeta) {
    // Figure out number of occupied orbitals
    int mult = mol.getMultiplicity();   //multiplicity
    int Ne = mol.getNElectrons();       //total electrons
    int Nd = Ne - (mult - 1);           //doubly occupied electrons
    if (Nd%2 != 0) MSG_FATAL("Invalid multiplicity");
    int Na = Nd/2 + (mult - 1);         //alpha electrons
    int Nb = Nd/2;                      //beta electrons

    // Make Fock operator contributions
    mrcpp::PoissonOperator P(*MRA, prec);
    mrcpp::ABGVOperator<3> D(*MRA, 0.5, 0.5);
    mrdft::XCFunctional xcfun(*MRA, false);
    xcfun.setFunctional("SLATERX");
    xcfun.setFunctional("VWN5C");
    xcfun.evalSetup(1);
    KineticOperator T(D);
    NuclearOperator V(mol.getNuclei(), prec);
    CoulombOperator J(&P);
    XCOperator XC(&xcfun);

    // Compute SAD density
    DensityVector rho_atomic;
    initial_guess::sad::project_atomic_densities(prec, mol, rho_atomic);
    Density &rho_xc = XC.getDensity(DENSITY::Total);
    Density &rho_j = J.getDensity();
    mrcpp::add(prec, rho_xc, rho_atomic);
    mrcpp::add(prec, rho_j, rho_atomic);

    // Project AO basis of hydrogen functions
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);

    // Compute Fock matrix
    Timer t_diag;
    Printer::printHeader(0, "Diagonalize AO Fock matrix");
    Timer t1;
    T.setup(prec);
    V.setup(prec);
    J.setup(prec);
    XC.setup(prec);
    ComplexMatrix t = T(Phi, Phi);
    ComplexMatrix v = V(Phi, Phi);
    ComplexMatrix j = J(Phi, Phi);
    ComplexMatrix xc = XC(Phi, Phi);
    ComplexMatrix F = S_m12.adjoint()*(t + v + j + xc)*S_m12;
    XC.clear();
    J.clear();
    V.clear();
    T.clear();
    t1.stop();
    Printer::printDouble(0, "Compute Fock matrix", t1.getWallTime(), 5);

    // Diagonalize Fock matrix
    Timer t2;
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(F.cols());
    es.compute(F);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    t2.stop();
    Printer::printDouble(0, "Diagonalize Fock matrix", t2.getWallTime(), 5);

    // Rotate orbitals and fill electrons by Aufbau
    Timer t3;
    OrbitalVector Psi;
    if (restricted) {
        for (int i = 0; i < Nd; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i = orbital::multiply(v_i, Phi, prec);
            psi_i.setSpin(SPIN::Paired);
            psi_i.setOcc(2);
            Psi.push_back(psi_i);
        }
        orbital::free(Phi);
    } else {
        OrbitalVector Psi_a;
        for (int i = 0; i < Na; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i = orbital::multiply(v_i, Phi, prec);
            psi_i.setSpin(SPIN::Alpha);
            psi_i.setOcc(1);
            Psi_a.push_back(psi_i);
        }
        OrbitalVector Psi_b;
        for (int i = 0; i < Nb; i++) {
            ComplexVector v_i = U.row(i);
            Orbital psi_i = orbital::multiply(v_i, Phi, prec);
            psi_i.setSpin(SPIN::Beta);
            psi_i.setOcc(1);
            Psi_b.push_back(psi_i);
        }
        orbital::free(Phi);
        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    t3.stop();
    Printer::printDouble(0, "Rotate orbitals", t3.getWallTime(), 5);

    t_diag.stop();
    Printer::printFooter(0, t_diag, 1);

    math_utils::print_matrix(0, es.eigenvalues(), "Eigenvalues", 10);

    return Psi;
}

void initial_guess::sad::project_atomic_densities(double prec,
                                                      const Molecule &mol,
                                                      DensityVector &rho_atomic) {
    Timer timer;
    Printer::printHeader(0, "Projecting Gaussian-type density");
    println(0, " Nr  Element         Rho_a         Rho_b         Rho_tot");
    Printer::printSeparator(0, '-');

    int oldprec = Printer::setPrecision(5);
    const Nuclei &nucs = mol.getNuclei();
    for (int k = 0; k < nucs.size(); k++) {
        const std::string& sym = nucs[k].getElement().getSymbol();

        std::stringstream bas;
        std::stringstream densa;
        std::stringstream densb;
        bas << "initial_guess/" << sym << ".bas";
        densa << "initial_guess/" << sym << ".densa";
        densb << "initial_guess/" << sym << ".densb";

        Density *rho_a = initial_guess::gto::project_density(prec, nucs[k], bas.str(), densa.str());
        Density *rho_b = initial_guess::gto::project_density(prec, nucs[k], bas.str(), densb.str());

        Density *rho_t = new Density(*MRA);
        mrcpp::copy_grid(*rho_t, *rho_a);
        mrcpp::copy_grid(*rho_t, *rho_b);
        mrcpp::add(-1.0, *rho_t, 1.0, *rho_a, 1.0, *rho_b);

        printout(0, std::setw(3) << k);
        printout(0, std::setw(7) << sym);
        printout(0, std::setw(22) << rho_a->integrate());
        printout(0, std::setw(14) << rho_b->integrate());
        printout(0, std::setw(14) << rho_t->integrate() << "\n");

        rho_atomic.push_back(std::make_tuple(1.0, rho_t));
    }
    Printer::setPrecision(oldprec);

    timer.stop();
    Printer::printFooter(0, timer, 2);
}


} //namespace mrchem

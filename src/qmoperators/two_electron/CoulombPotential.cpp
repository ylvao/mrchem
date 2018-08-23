#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombPotential.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/density_utils.h"

using mrcpp::FunctionTree;
using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */
CoulombPotential::CoulombPotential(PoissonOperator *P, OrbitalVector *Phi,
                                   OrbitalVector *Phi_1, int order)
        : QMPotential(1), density(*MRA), density_1(*MRA), orbitals(Phi),
          orbitals_1(Phi), poisson(P) {
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application of the Poisson
 * operator to the density. If the density is not available it is computed
 * from the current orbitals (assuming that the orbitals are available).
 * For first-order perturbations the first order density and the Hessian will be
 * computed. In order to make the Hessian available to CoulombOperator, it is stored in the 
 * potential function instead of the zeroth-order potential.
 *
 */
void CoulombPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity(prec);
    if(this->order =  1) {
        setupPotential(prec);
    }
    else if (this->order = 2) {
        setupDensity_1(prec);
        setupHessian(prec);
    }
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    QMFunction::free(); // delete FunctionTree pointers
    clearApplyPrec();   // apply_prec = -1
    mrcpp::clear_grid(this->density.real()); // clear MW coefs but keep the grid
}

/** @brief compute electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the electron density as the sum of squares of the orbitals.
 */
void CoulombPotential::setupDensity(double prec) {
    if (hasDensity()) return;
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    Density &rho = this->density;

    Timer timer;
    density::compute(prec, rho, Phi, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

/** @brief compute first-order electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the first-order perturbed electron density.
 */
void CoulombPotential::setupPertDensity(double prec) {
    if (hasDensity_1()) return;
    if (not hasDensity()) MSG_ERROR("Ground-state density not initialized");
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_1 == nullptr) MSG_ERROR("Perturbed orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Phi_1 = *this->orbitals_1;
    Density &rho = this->density;
    Density &rho_1 = this->density_1;

    Timer timer;
    density::compute(prec, rho_1, Phi, Phi_1, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

/** @brief compute first-order electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the first-order perturbed electron density.
 */
void CoulombPotential::setupDensity(double prec) {
    if (hasDensity_1()) return;
    if (not hasDensity()) MSG_ERROR("Ground-state density not initialized");
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->orbitals_1 == nullptr) MSG_ERROR("Perturbed orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &Phi_1 = *this->orbitals;
    Density &rho = this->density;
    Density &rho_1 = this->density_1;

    Timer timer;
    density::compute(prec, rho, Phi, DENSITY::Total);
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Coulomb density", n, t);
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
void CoulombPotential::setupPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");
    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    PoissonOperator &P = *this->poisson;
    QMPotential &V = *this;
    Density &rho = this->density;

    // Adjust precision by system size
    double abs_prec = prec/rho.real().integrate();

    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::apply(abs_prec, V.real(), P, rho.real());
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);
}

} //namespace mrchem

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "CoulombHessian.h"
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
CoulombHessian::CoulombHessian(OrbitalVector *X) {
    this->pertX = X;
    this->pertY = nullptr;
    this->density.alloc(NUMBER::Real);
}

CoulombHessian::CoulombHessian(OrbitalVector *X,
                               OrbitalVector *Y) {
    this->pertX = X;
    this->pertY = Y;
    this->density.alloc(NUMBER::Real);
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
void CoulombHessian::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    setupDensity(prec);
    setupHessian(prec);
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombHessian::clear() {
    QMFunction::free(); // delete FunctionTree pointers
    clearApplyPrec();   // apply_prec = -1
    mrcpp::clear_grid(this->density.real()); // clear MW coefs but keep the grid
}

/** @brief compute first-order electron density
 *
 * @param[in] prec: apply precision
 *
 * This will compute the first-order perturbed electron density.
 */
void CoulombHessian::setupDensity(double prec) {
    if (not hasDensity()) MSG_ERROR("Ground-state density not initialized");
    if (this->orbitals == nullptr) MSG_ERROR("Orbitals not initialized");
    if (this->pertX == nullptr) MSG_ERROR("Perturbed orbitals not initialized");

    OrbitalVector &Phi = *this->orbitals;
    OrbitalVector &X = *this->pertX;
    OrbitalVector &Y = *this->pertY;
    Density &rho = this->density;

    Timer timer;
    if(this->pertY == nullptr) { // For static perturbations pertY is not required/defined
        density::compute(prec, rho, Phi, X, DENSITY::Total);
    } else {
        density::compute(prec, rho, Phi, X, Y, DENSITY::Total);
    }
    timer.stop();
    double t = timer.getWallTime();
    int n = rho.getNNodes();
    Printer::printTree(0, "Perturbed Coulomb density", n, t);
}

/** @brief compute Coulomb hessian
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb hessian by application o the Poisson operator
 * to the precomputed, perturbed electron density.
 */
void CoulombHessian::setupHessian(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");
    if (hasReal()) MSG_ERROR("Hessian not properly cleared");
    if (hasImag()) MSG_ERROR("Hessian not properly cleared");

    PoissonOperator &P = *this->poisson;
    QMPotential &V = *this;
    Density &rho = this->density;

    /*    
    // Adjust precision by system size
    double abs_prec = prec/rho.real().integrate();
    */
    
//LUCA: this is currently using only the real part of the perturbed density. Should be OK for now
    Timer timer;
    V.alloc(NUMBER::Real);
    mrcpp::apply(prec, V.real(), P, rho.real());
    timer.stop();
    int n = V.getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb hessian", n, t);
}

} //namespace mrchem

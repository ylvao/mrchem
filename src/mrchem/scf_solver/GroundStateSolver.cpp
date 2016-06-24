#include <Eigen/Eigenvalues>

#include "GroundStateSolver.h"
#include "HelmholtzOperatorSet.h"
#include "FockOperator.h"
#include "PoissonOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombPotential.h"
#include "ExchangePotential.h"
#include "XCPotential.h"
#include "XCFunctional.h"
#include "SCFEnergy.h"
//#include "KAIN.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

GroundStateSolver::GroundStateSolver(const MultiResolutionAnalysis<3> &mra,
                                     HelmholtzOperatorSet &h)
        : SCF(mra, h),
          fOper_n(0),
          fMat_n(0),
          orbitals_n(0),
          orbitals_np1(0),
          dOrbitals_n(0) {
}

GroundStateSolver::~GroundStateSolver() {
    if (this->fOper_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->fMat_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->orbitals_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->orbitals_np1 != 0) MSG_ERROR("Solver not properly cleared");
    if (this->dOrbitals_n != 0) MSG_ERROR("Solver not properly cleared");
}

/** Computes the Helmholtz argument for the i-th orbital.
 *
 * Argument contains the potential operator acting on orbital i, and the sum
 * of all orbitals weighted by the Fock matrix. The effect of using inexact
 * Helmholtz operators are included in Lambda, wich is a diagonal matrix
 * with the actual lambda parameters used in the Helmholtz operators.
 *
 * greenArg = \hat{V}orb_i + \sum_j (\Lambda_{ij}-F_{ij})orb_j
 */
Orbital* GroundStateSolver::getHelmholtzArgument(int i,
                                                 MatrixXd &F,
                                                 OrbitalVector &phi,
                                                 bool adjoint) {
    FockOperator &fock = *this->fOper_n;

    double coef = -1.0/(2.0*pi);
    Orbital &phi_i = phi.getOrbital(i);

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd LmF = L - F;

    Orbital *part_1 = fock.applyPotential(phi_i);
    Orbital *part_2 = calcMatrixPart(i, LmF, phi);

    if (part_1 == 0) part_1 = new Orbital(phi_i);
    if (part_2 == 0) part_2 = new Orbital(phi_i);

    Timer timer;
    timer.restart();
    Orbital *arg = new Orbital(phi_i);
    this->add(*arg, coef, *part_1, coef, *part_2);

    double time = timer.getWallTime();
    int nNodes = arg->getNNodes();
    TelePrompter::printTree(2, "Added arguments", nNodes, time);

    if (part_1 != 0) delete part_1;
    if (part_2 != 0) delete part_2;

    return arg;
}

double GroundStateSolver::calcProperty() {
    MatrixXd &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi = *this->orbitals_n;

    SCFEnergy scfEnergy;
    scfEnergy.compute(fock, F, phi);
    this->energy.push_back(scfEnergy);
    return scfEnergy.getElectronicEnergy();
}

double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return fabs(getUpdate(this->property, iter, false));
}

void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();

    TelePrompter::printHeader(0, "                    Energy                 Update      Done ");
    printUpdate(" Kinetic  ",  T_1,  T_1 -  T_0);
    printUpdate(" N-E      ",  V_1,  V_1 -  V_0);
    printUpdate(" Coulomb  ",  J_1,  J_1 -  J_0);
    printUpdate(" Exchange ",  K_1,  K_1 -  K_0);
    printUpdate(" X-C      ", XC_1, XC_1 - XC_0);
    TelePrompter::printSeparator(0, '-');
    printUpdate(" Total    ",  E_1,  E_1 -  E_0);
    TelePrompter::printSeparator(0, '=');
}

/** Minimize the spatial extension of orbitals, by a transformation of orbitals
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 *	which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 */
void GroundStateSolver::localize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    NOT_IMPLEMENTED_ABORT;
//    double optTime = 0.0;
//    double rotTime = 0.0;
//    double totTime = 0.0;

//    mpi::timer totTimer, timer;
//    totTimer.restart();

//    int oldPrec = TelePrompter::setPrecision(5);
//    println(0, endl);

//    timer.restart();
//    RR rr(*this);
//    int nIter = rr.maximize();//compute total U, rotation matrix
//    optTime += timer.elapsed();
//    println(0, "Localized after iteration               " << setw(20) << nIter);

//    MatrixXd U;
//    if (nIter > 0) {
//        U = rr.getTotalU().transpose();
//    } else {
//        timer.restart();
//        U = calcOrthonormalizationMatrix();
//        optTime += timer.elapsed();
//    }

//    println(0, "Calculating rotation matrix                      " << optTime);
//    println(0, "Total time localization                          " << totTimer.elapsed());
//    printout(0, endl);

//    TelePrompter::setPrecision(oldPrec);
//    return U;
}

/** Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * This operation includes the orthonormalization using the overlap matrix.
 */
void GroundStateSolver::diagonalize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    MatrixXd S_m12 = calcOrthonormalizationMatrix(phi);
    F = S_m12.transpose()*F*S_m12;

    Timer timer;
    timer.restart();
    printout(1, "Calculating diagonalization matrix               ");

    SelfAdjointEigenSolver<MatrixXd> es(F.cols());
    es.compute(F);
    MatrixXd M = es.eigenvectors();
    MatrixXd U = M.transpose()*S_m12;

    println(1, timer.getWallTime());

    F = es.eigenvalues().asDiagonal();
    fock.rotate(U);
    this->add.rotate(phi, U);
}

void GroundStateSolver::orthonormalize(FockOperator &fock, MatrixXd &F, OrbitalVector &phi) {
    MatrixXd U = calcOrthonormalizationMatrix(phi);
    F = U.transpose()*F*U;
    fock.rotate(U);
    this->add.rotate(phi, U);
}

MatrixXd GroundStateSolver::calcOrthonormalizationMatrix(OrbitalVector &phi) {
    Timer timer;
    timer.restart();
    printout(1, "Calculating orthonormalization matrix            ");

    MatrixXd S_tilde = phi.calcOverlapMatrix().real();
    SelfAdjointEigenSolver<MatrixXd> es(S_tilde.cols());
    es.compute(S_tilde);

    MatrixXd A = es.eigenvalues().asDiagonal();
    for (int i = 0; i < A.cols(); i++) {
        A(i,i) = pow(A(i,i), -1.0/2.0);
    }
    MatrixXd B = es.eigenvectors();
    MatrixXd U = B*A*B.transpose();

    println(1, timer.getWallTime());
    return U;
}

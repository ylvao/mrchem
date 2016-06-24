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
                                     HelmholtzOperatorSet &h,
                                     KAIN *k)
        : SCF(mra, h),
          fOper_n(0),
          fOper_np1(0),
          orbitals_n(0),
          orbitals_np1(0),
          dOrbitals_n(0),
          fMat_n(0),
          fMat_np1(0),
          dfMat_n(0),
          kain(k) {
}

GroundStateSolver::~GroundStateSolver() {
    if (this->orbitals_np1 != 0) MSG_ERROR("Solver not properly cleared");
    if (this->dOrbitals_n != 0) MSG_ERROR("Solver not properly cleared");
    if (this->fMat_np1 != 0) MSG_ERROR("Solver not properly cleared");
    if (this->dfMat_n != 0) MSG_ERROR("Solver not properly cleared");
    this->kain = 0;
}

void GroundStateSolver::setup(FockOperator &f_oper,
                              MatrixXd &f_mat,
                              OrbitalVector &phi) {
    this->orbitals_n = &phi;
    this->fMat_n = &f_mat;
    this->fOper_n = &f_oper;

    this->orbitals_np1 = new OrbitalVector(*this->orbitals_n);
    this->dOrbitals_n = new OrbitalVector(*this->orbitals_n);

    this->fMat_np1 = new MatrixXd;
    this->dfMat_n = new MatrixXd;
}

void GroundStateSolver::clear() {
    clearUpdates();
    this->nIter = 0;
    this->orbitals_n = 0;
    this->fMat_n = 0;
    this->fOper_n = 0;

    if (this->orbitals_np1 != 0) delete this->orbitals_np1;
    if (this->dOrbitals_n != 0) delete this->dOrbitals_n;

    if (this->fMat_np1 != 0) delete this->fMat_np1;
    if (this->dfMat_n != 0) delete this->dfMat_n;

    this->fMat_np1 = 0;
    this->dfMat_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

//    if (this->kain != 0) this->kain->clear();
    resetPrecision();
}

void GroundStateSolver::setupUpdates() {
    NOT_IMPLEMENTED_ABORT;
//    CoulombOperator *j_np1 = 0;
//    ExchangeOperator *k_np1 = 0;
//    XCOperator *xc_np1 = 0;

//    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
//    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
//    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
//    XCOperator *xc_n = this->fOper_n->getXCOperator();

//    if (j_n != 0) j_np1 = new CoulombPotential(*j_n, *this->phi_np1);
//    if (k_n != 0) k_np1 = new ExchangePotential(*k_n, *this->phi_np1);
//    if (xc_n != 0) xc_np1 = new XCPotential(*xc_n, *this->phi_np1);

//    this->fOper_np1 = new FockOperator(0, v_n, j_np1, k_np1, xc_np1);
}

void GroundStateSolver::clearUpdates() {
    if (this->fOper_np1 != 0) {
        CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
        ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
        XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

        if (j_np1 != 0) delete j_np1;
        if (k_np1 != 0) delete k_np1;
        if (xc_np1 != 0) delete xc_np1;

        delete this->fOper_np1;
    }
    this->fOper_np1 = 0;
}

bool GroundStateSolver::optimize() {
    int oldPrec = TelePrompter::setPrecision(5);
    bool converged = optimizeOrbitals();
    if (converged and getPropertyThreshold() > 0.0) {
        int oldRot = this->iterPerRotation;
        setRotation(0);
        setupUpdates();
        converged = optimizeEnergy();
        clearUpdates();
        setRotation(oldRot);
    }
    printConvergence(converged);
    TelePrompter::setPrecision(oldPrec);
    return converged;
}

bool GroundStateSolver::optimizeOrbitals() {
    MatrixXd &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;
    adjustPrecision(err_o);

    fock.setup(getOrbitalPrecision());
    F = fock(phi_n, phi_n);

    bool converged = false;
    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
        Timer timer;
        timer.restart();

        {   // Initialize SCF cycle
            printCycle();
            adjustPrecision(err_o);
        }

        {   // Localize/diagonalize/orthogonalize
            MatrixXd U;
            if (needLocalization()) {
                U = calcLocalizationMatrix(phi_n);
            } else if (needDiagonalization()) {
                U = calcDiagonalizationMatrix(phi_n, F);
//                if (this->kain != 0) this->kain->clear();
            } else {
                U = calcOrthonormalizationMatrix(phi_n);
            }
            fock.rotate(U);
            F = U.transpose()*F*U;
            this->add.rotate(phi_n, U);
        }

        {   // Compute electronic energy
            double prop = calcProperty();
            this->property.push_back(prop);
        }

        {   // Iterate Helmholtz operators
            this->helmholtz->initialize(F.diagonal());
            applyHelmholtzOperators(phi_np1, phi_n, F);
            fock.clear();
        }

        {   // Orthonormalize
            MatrixXd U = calcOrthonormalizationMatrix(phi_np1);
            this->add.rotate(phi_np1, U);
        }

        {   // Compute orbital updates
            this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n);
            phi_np1.clear();
        }

        { // Employ KAIN accelerator
//            if (this->kain != 0) this->kain->pushBack(phi_n, dPhi_n);
//            if (this->kain != 0) this->kain->calcUpdates(phi_n, dPhi_n);
        }

//        printTreeSizes();

        {   // Compute errors
            VectorXd errors = dPhi_n.getNorms();
            phi_n.setErrors(errors);
            err_o = errors.maxCoeff();
            err_t = sqrt(errors.dot(errors));
            this->orbError.push_back(err_t);
        }

        {   // Update orbitals
            this->add.inPlace(phi_n, 1.0, dPhi_n);
            dPhi_n.clear();
        }

        {   // Orthonormalize
            MatrixXd U = calcOrthonormalizationMatrix(phi_n);
            this->add.rotate(phi_n, U);
//            phi_n.orthogonalize();
//            phi_n.normalize();
        }

        {   // Compute Fock matrix
            fock.setup(getOrbitalPrecision());
            F = fock(phi_n, phi_n);
        }

        {   // Finalize SCF cycle
            printOrbitals(F, phi_n);
            printProperty();
            printTimer(timer.getWallTime());
        }

        if (err_o < getOrbitalThreshold()) {
            converged = true;
            break;
        }
    }
//    if (this->kain != 0) this->kain->clear();
    fock.clear();
    return converged;
}

bool GroundStateSolver::optimizeEnergy() {
    NOT_IMPLEMENTED_ABORT;
//    boost::timer scfTimer;
//    double err_p = 1.0;
//    double err_o = this->phi_n->getErrors().maxCoeff();
//    double err_t = 1.0;

//    bool first = true;
//    bool converged = false;
//    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
//        scfTimer.restart();
//        printCycle();
//        adjustPrecision(err_o);

//        printMatrix(1, this->phi_n->calcOverlapMatrix(), 'S');
//        printMatrix(1, *this->fMat_n, 'F');

//        this->fOper_n->setup(getOrbitalPrecision());
//        double prop = calcProperty();
//        this->property.push_back(prop);

//        this->helmholtz->initialize(this->fMat_n->diagonal());
//        applyHelmholtzOperators(*this->phi_np1, *this->fMat_n, *this->phi_n);
//        calcOrbitalUpdates();
//        calcFockMatrixUpdate();
//        printTreeSizes();
//        this->dPhi_n->clear();

//        *this->fMat_np1 = *this->fMat_n + *this->dfMat_n;

//        printMatrix(1, *this->fMat_np1, 'F');
//        rotate(*this->phi_np1, *this->fMat_np1);
//        printMatrix(1, *this->fMat_np1, 'F');

//        calcOrbitalUpdates();
//        *this->dfMat_n = *this->fMat_np1 - *this->fMat_n;
//        this->phi_np1->clear();

//        err_p = calcPropertyError();
//        err_o = calcOrbitalError();
//        err_t = calcTotalError();
//        this->orbError.push_back(err_t);

//        this->phi_n->addInPlace(*this->dPhi_n);
//        *this->fMat_n += *this->dfMat_n;

//        this->phi_n->orthonormalize(getOrbitalPrecision(), this->fMat_n);
//        this->phi_n->setErrors(this->dPhi_n->getNorms());

//        this->fOper_n->clear();
//        this->dPhi_n->clear();

//        printOrbitals(*this->fMat_n, *this->phi_n);
//        printProperty();
//        printTimer(scfTimer.elapsed());

//        if (err_p < getPropertyThreshold() and not first) {
//            converged = true;
//            break;
//        }
//        first = false;
//    }
//    return converged;
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

double GroundStateSolver::calcTotalError() const {
    VectorXd errors = this->dOrbitals_n->getNorms();
    return sqrt(errors.dot(errors));
}

double GroundStateSolver::calcOrbitalError() const {
    return this->dOrbitals_n->getNorms().maxCoeff();
}

double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return fabs(getUpdate(this->property, iter, false));
}

//void GroundStateSolver::rotate(OrbitalVector &phi,
//                               MatrixXd &f_mat,
//                               FockOperator *f_oper) {
//        NOT_IMPLEMENTED_ABORT;
//    Timer timer;
//    timer.restart();
//    TelePrompter::printHeader(1,"Rotating Orbitals");

//    MatrixXd U;
//    if (needLocalization()) {
//        U = phi.localize(getOrbitalPrecision(), &f_mat);
//    } else if (needDiagonalization()) {
//        U = phi.diagonalize(getOrbitalPrecision(), &f_mat);
//        if (this->acc != 0) this->acc->clear();
//    } else {
//        U = phi.orthonormalize(getOrbitalPrecision(), &f_mat);
//    }

//    if (f_oper != 0) {
//        f_oper->rotate(getOrbitalPrecision(), U);
//    }
//    TelePrompter::printFooter(1, timer, 2);
//}

void GroundStateSolver::calcOrbitalUpdates() {
    NOT_IMPLEMENTED_ABORT;
//    boost::timer timer;
//    timer.restart();
//    println(1,"                                                            ");
//    println(1,"-------------------- Calculating Updates -------------------");
//    this->dPhi_n->clear();
//    this->dPhi_n->add(1.0, *this->phi_np1, -1.0, *this->phi_n);
//    printout(1,"---------------- Elapsed time: " << timer.elapsed());
//    println(1," -----------------");
//    println(1,"                                                            ");
}

void GroundStateSolver::calcFockMatrixUpdate() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");
//    boost::timer timer, totTimer;
//    println(0,"                                                            ");
//    println(0,"=============== Computing Fock Matrix update ===============");
//    println(0,"                                                            ");

//    timer.restart();
//    MatrixXd dS_1 = this->dPhi_n->calcOverlapMatrix(*this->phi_n);
//    MatrixXd dS_2 = this->phi_np1->calcOverlapMatrix(*this->dPhi_n);
//    double time_s = timer.elapsed();
//    TelePrompter::printDouble(0, "Overlap matrices", time_s);

//    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
//    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
//    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
//    XCOperator *xc_n = this->fOper_n->getXCOperator();

//    // Nuclear potential matrix is computed explicitly
//    timer.restart();
//    MatrixXd dV_n = (*v_n)(*this->phi_np1, *this->dPhi_n);
//    double time_v = timer.elapsed();
//    TelePrompter::printDouble(0, "Nuclear potential matrix", time_v);

//    timer.restart();
//    FockOperator f_n(0, 0, j_n, k_n, xc_n);
//    MatrixXd F_n = f_n(*this->phi_np1, *this->phi_n);
//    double time_f_n = timer.elapsed();
//    TelePrompter::printDouble(0, "Fock matrix n", time_f_n);
//    this->fOper_n->clear();

//    // The n+1 Fock operator needs orthonormalized orbitals
//    this->phi_np1->orthonormalize(getOrbitalPrecision());

//    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
//    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
//    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

//    println(0,"                                                            ");
//    // Do not setup exchange, it must be applied on the fly anyway
//    if (j_np1 != 0) j_np1->setup(getOrbitalPrecision());
//    if (k_np1 != 0) j_np1->QMOperator::setup(getOrbitalPrecision());
//    if (xc_np1 != 0) xc_np1->setup(getOrbitalPrecision());
//    println(0,"                                                            ");

//    // Computing potential matrix excluding nuclear part
//    timer.restart();
//    FockOperator f_np1(0, 0, j_np1, k_np1, xc_np1);
//    MatrixXd F_1 = f_np1(*this->phi_n, *this->phi_n);
//    MatrixXd F_2 = f_np1(*this->phi_n, *this->dPhi_n);
//    MatrixXd F_np1 = F_1 + F_2 + F_2.transpose();
//    //MatrixXd F_3 = f_np1(*this->dPhi_n, *this->phi_n);
//    //MatrixXd F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
//    //MatrixXd F_np1 = F_1 + F_2 + F_3 + F_4;
//    double time_f_np1 = timer.elapsed();
//    TelePrompter::printDouble(0, "Fock matrix n+1", time_f_np1);

//    f_np1.clear();
//    this->phi_np1->clear();
//    this->phi_np1->add(1.0, *this->phi_n, 1.0, *this->dPhi_n);

//    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
//    MatrixXd dF_n = F_np1 - F_n;
//    MatrixXd dF_1 = dS_1*(*this->fMat_n);
//    MatrixXd dF_2 = dS_2*L;

//    // Adding up the pieces
//    *this->dfMat_n = dF_1 + dF_2 + dV_n + dF_n;

//    //Symmetrizing
//    MatrixXd sym = *this->dfMat_n + this->dfMat_n->transpose();
//    *this->dfMat_n = 0.5 * sym;

//    double t = totTimer.elapsed();
//    println(0, "                                                            ");
//    println(0, "================ Elapsed time: " << t << " =================");
//    println(0, "                                                            ");
//    println(0, "                                                            ");
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
                                                 OrbitalVector &phi,
                                                 MatrixXd &F,
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

/** Prints the number of trees and nodes kept in the solver at the given moment */
int GroundStateSolver::printTreeSizes() const {
    NOT_IMPLEMENTED_ABORT;
//    println(0,"                                                            ");
//    println(0,"                                                            ");
//    println(0,"------------------- Printing Tree sizes --------------------");
//    println(0,"                                                            ");

//    int n2D = 0;
//    int n3D = 0;
//    if (this->helmholtz != 0) n2D += this->helmholtz->printTreeSizes();
//    if (this->fOper_n != 0) n3D += this->fOper_n->printTreeSizes();
//    if (this->fOper_np1 != 0) n3D += this->fOper_np1->printTreeSizes();
//    if (this->phi_n != 0) n3D += this->phi_n->printTreeSizes();
//    if (this->phi_np1 != 0) n3D += this->phi_np1->printTreeSizes();
//    if (this->dPhi_n != 0) n3D += this->dPhi_n->printTreeSizes();
//    if (this->acc != 0) n3D += this->acc->printTreeSizes();

//    println(0,"                                                            ");
//    println(0," Total number of 2D nodes                " << setw(18) << n2D);
//    println(0," Total number of 3D nodes                " << setw(18) << n3D);
//    println(0,"                                                            ");
//    println(0,"------------------------------------------------------------");
//    println(0,"                                                            ");

//    return n3D;
}

/** Minimize the spatial extension of orbitals, by a transformation of orbitals
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 *	which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * For details see the tex documentation in doc directory
 */
MatrixXd GroundStateSolver::calcLocalizationMatrix(OrbitalVector &phi) {
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
MatrixXd GroundStateSolver::calcDiagonalizationMatrix(OrbitalVector &phi, MatrixXd &F_tilde) {
    Timer timer;
    timer.restart();
    printout(1, "Calculating diagonalization matrix               ");

    MatrixXd S_m12 = calcOrthonormalizationMatrix(phi);
    MatrixXd F = S_m12*F_tilde*S_m12.transpose();

    SelfAdjointEigenSolver<MatrixXd> es(F.cols());
    es.compute(F);
    MatrixXd M = es.eigenvalues();
    MatrixXd U = M.transpose()*S_m12;

    println(1, timer.getWallTime());
    return U;
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

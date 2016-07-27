#include "EnergyOptimizer.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "HelmholtzOperatorSet.h"
#include "eigen_disable_warnings.h"

using namespace std;
using namespace Eigen;

EnergyOptimizer::EnergyOptimizer(const MultiResolutionAnalysis<3> &mra,
                                 HelmholtzOperatorSet &h)
        : GroundStateSolver(mra, h),
          fOper_np1(0) {
}

EnergyOptimizer::~EnergyOptimizer() {
    if (this->fOper_np1 != 0) MSG_ERROR("Solver not properly cleared");
}

void EnergyOptimizer::setup(FockOperator &fock, OrbitalVector &phi, MatrixXd &F,
                            FockOperator &fock_np1, OrbitalVector &phi_np1) {
    this->fMat_n = &F;
    this->fOper_n = &fock;
    this->fOper_np1 = &fock_np1;

    this->orbitals_n = &phi;
    this->orbitals_np1 = &phi_np1;
    this->dOrbitals_n = new OrbitalVector(phi);
}

void EnergyOptimizer::clear() {
    this->nIter = 0;
    this->fMat_n = 0;
    this->fOper_n = 0;
    this->fOper_np1 = 0;

    delete this->dOrbitals_n;

    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

    resetPrecision();
}

bool EnergyOptimizer::optimize() {
    MatrixXd &F_n = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;
    double err_p = 1.0;

    bool converged = false;
    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
        Timer timer;
        timer.restart();

        // Initialize SCF cycle
        printCycle();
        adjustPrecision(err_o);

        // Compute electronic energy
        fock.setup(getOrbitalPrecision());
        double E = calcProperty();
        this->property.push_back(E);

        // Iterate Helmholtz operators
        this->helmholtz->initialize(F_n.diagonal());
        applyHelmholtzOperators(phi_np1, F_n, phi_n);
        this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n);

        // Compute errors
        VectorXd errors = dPhi_n.getNorms();
        phi_n.setErrors(errors);
        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
        err_p = calcPropertyError();
        this->orbError.push_back(err_t);

        // Compute Fock matrix
        MatrixXd F_np1 = F_n + calcFockMatrixUpdate();
        phi_n.clear();
        dPhi_n.clear();
        fock.clear();

        // Rotate orbitals
        if (needLocalization()) {
            localize(fock, F_np1, phi_np1);
        } else if (needDiagonalization()) {
            diagonalize(fock, F_np1, phi_np1);
        } else {
            orthonormalize(fock, F_np1, phi_np1);
        }

        // Update orbitals and Fock matrix
        int nOrbs = phi_np1.size();
        MatrixXd U = MatrixXd::Identity(nOrbs, nOrbs);
        F_n = U.transpose()*F_np1*U;
        this->add.rotate(phi_n, U, phi_np1);
        phi_np1.clear();

        // Finalize SCF cycle
        printOrbitals(F_n, phi_n);
        printProperty();
        printTimer(timer.getWallTime());

        if (err_p < getPropertyThreshold()) {
            converged = true;
            break;
        }
    }
    printConvergence(converged);
    return converged;
}

MatrixXd EnergyOptimizer::calcFockMatrixUpdate() {
    if (this->fOper_np1 == 0) MSG_FATAL("Operator not initialized");

    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    TelePrompter::printHeader(0,"Computing Fock matrix update");

    Timer timer;
    timer.restart();

    MatrixXd dS_1 = dPhi_n.calcOverlapMatrix(phi_n).real();
    MatrixXd dS_2 = phi_np1.calcOverlapMatrix(dPhi_n).real();

    NuclearPotential *v_n = this->fOper_n->getNuclearPotential();
    CoulombOperator *j_n = this->fOper_n->getCoulombOperator();
    ExchangeOperator *k_n = this->fOper_n->getExchangeOperator();
    XCOperator *xc_n = this->fOper_n->getXCOperator();

    MatrixXd dV_n;
    {   // Nuclear potential matrix is computed explicitly
        Timer timer;
        timer.restart();
        dV_n = (*v_n)(phi_np1, dPhi_n);
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Nuclear potential matrix", t);
    }

    MatrixXd F_n;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        timer.restart();
        FockOperator fock_n(this->MRA, 0, 0, j_n, k_n, xc_n);
        F_n = fock_n(phi_np1, phi_n);
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n", t);
    }

    {   // The n+1 Fock operator needs orthonormalized orbitals
        MatrixXd U = calcOrthonormalizationMatrix(phi_np1);
        this->add.rotate(phi_np1, U);
    }

    CoulombOperator *j_np1 = this->fOper_np1->getCoulombOperator();
    ExchangeOperator *k_np1 = this->fOper_np1->getExchangeOperator();
    XCOperator *xc_np1 = this->fOper_np1->getXCOperator();

    println(0,"                                                            ");
    // Do not setup exchange, it must be applied on the fly anyway
    if (j_np1 != 0) j_np1->setup(getOrbitalPrecision());
    if (k_np1 != 0) k_np1->QMOperator::setup(getOrbitalPrecision());
    if (xc_np1 != 0) xc_np1->setup(getOrbitalPrecision());
    println(0,"                                                            ");

    MatrixXd F_np1;
    {   // Computing potential matrix excluding nuclear part
        Timer timer;
        timer.restart();

        FockOperator fock_np1(this->MRA, 0, 0, j_np1, k_np1, xc_np1);
        MatrixXd F_1 = fock_np1(phi_n, phi_n);
        MatrixXd F_2 = fock_np1(phi_n, dPhi_n);
        fock_np1.clear();

        F_np1 = F_1 + F_2 + F_2.transpose();
        //MatrixXd F_3 = f_np1(*this->dPhi_n, *this->phi_n);
        //MatrixXd F_4 = f_np1(*this->dPhi_n, *this->dPhi_n);
        //MatrixXd F_np1 = F_1 + F_2 + F_3 + F_4;
        double t = timer.getWallTime();
        TelePrompter::printDouble(0, "Fock matrix n+1", t);
    }

    // Re-computing non-orthogonal phi_np1
    phi_np1.clear();
    this->add(phi_np1, 1.0, phi_n, 1.0, dPhi_n);

    MatrixXd L = this->helmholtz->getLambda().asDiagonal();
    MatrixXd dF_1 = dS_1*(*this->fMat_n);
    MatrixXd dF_2 = dS_2*L;
    MatrixXd dF_3 = F_np1 - F_n;

    // Adding up the pieces
    MatrixXd dF_n = dV_n + dF_1 + dF_2 + dF_3;

    // Symmetrizing
    MatrixXd sym = dF_n + dF_n.transpose();
    dF_n = 0.5 * sym;

    TelePrompter::printFooter(0, timer, 2);
    return dF_n;
}

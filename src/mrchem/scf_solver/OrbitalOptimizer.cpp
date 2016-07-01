#include "OrbitalOptimizer.h"
#include "OrbitalVector.h"
#include "FockOperator.h"
#include "HelmholtzOperatorSet.h"
#include "Accelerator.h"

using namespace std;
using namespace Eigen;

OrbitalOptimizer::OrbitalOptimizer(const MultiResolutionAnalysis<3> &mra,
                                   HelmholtzOperatorSet &h,
                                   Accelerator *k)
        : GroundStateSolver(mra, h),
          kain(k) {
}

OrbitalOptimizer::~OrbitalOptimizer() {
    this->kain = 0;
}

void OrbitalOptimizer::setup(FockOperator &fock,
                              OrbitalVector &phi,
                              MatrixXd &F) {
    this->fMat_n = &F;
    this->fOper_n = &fock;

    this->orbitals_n = &phi;
    this->orbitals_np1 = new OrbitalVector(phi);
    this->dOrbitals_n = new OrbitalVector(phi);
}

void OrbitalOptimizer::clear() {
    this->nIter = 0;
    this->fMat_n = 0;
    this->fOper_n = 0;

    delete this->orbitals_np1;
    delete this->dOrbitals_n;

    this->orbitals_n = 0;
    this->orbitals_np1 = 0;
    this->dOrbitals_n = 0;

    if (this->kain != 0) this->kain->clear();
    resetPrecision();
}

bool OrbitalOptimizer::optimize() {
    MatrixXd &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &phi_n = *this->orbitals_n;
    OrbitalVector &phi_np1 = *this->orbitals_np1;
    OrbitalVector &dPhi_n = *this->dOrbitals_n;

    double err_o = phi_n.getErrors().maxCoeff();
    double err_t = 1.0;

    fock.setup(getOrbitalPrecision());
    F = fock(phi_n, phi_n);

    bool converged = false;
    while(this->nIter++ < this->maxIter or this->maxIter < 0) {
        Timer timer;
        timer.restart();

        // Initialize SCF cycle
        printCycle();
        adjustPrecision(err_o);

        // Rotate orbitals
        if (needLocalization()) {
            localize(fock, F, phi_n);
        } else if (needDiagonalization()) {
            diagonalize(fock, F, phi_n);
            if (this->kain != 0) this->kain->clear();
        } else {
            orthonormalize(fock, F, phi_n);
        }

        // Compute electronic energy
        double E = calcProperty();
        this->property.push_back(E);

        // Iterate Helmholtz operators
        this->helmholtz->initialize(F.diagonal());
        applyHelmholtzOperators(phi_np1, F, phi_n);
        fock.clear();

        orthonormalize(fock, F, phi_np1);

        // Compute orbital updates
        this->add(dPhi_n, 1.0, phi_np1, -1.0, phi_n);
        phi_np1.clear();

        // Employ KAIN accelerator
        if (this->kain != 0) {
            this->kain->setPrecision(this->orbPrec[0]);
            this->kain->accelerate(phi_n, dPhi_n);
        }

        // Compute errors
        VectorXd errors = dPhi_n.getNorms();
        phi_n.setErrors(errors);
        err_o = errors.maxCoeff();
        err_t = sqrt(errors.dot(errors));
        this->orbError.push_back(err_t);

        // Update orbitals
        this->add.inPlace(phi_n, 1.0, dPhi_n);
        dPhi_n.clear();

        orthonormalize(fock, F, phi_n);

        // Compute Fock matrix
        fock.setup(getOrbitalPrecision());
        F = fock(phi_n, phi_n);

        // Finalize SCF cycle
        printOrbitals(F, phi_n);
        printProperty();
        printTimer(timer.getWallTime());

        if (err_o < getOrbitalThreshold()) {
            converged = true;
            break;
        }
    }
    if (this->kain != 0) this->kain->clear();
    fock.clear();
    printConvergence(converged);
    return converged;
}

/** Prints the number of trees and nodes kept in the solver at the given moment */
void OrbitalOptimizer::printTreeSizes() const {
    TelePrompter::printHeader(0, "Printing Tree Sizes");

    int nNodes = 0;
    if (this->fOper_n != 0) nNodes += this->fOper_n->printTreeSizes();
    if (this->orbitals_n != 0) nNodes += this->orbitals_n->printTreeSizes();
    if (this->orbitals_np1 != 0) nNodes += this->orbitals_np1->printTreeSizes();
    if (this->dOrbitals_n != 0) nNodes += this->dOrbitals_n->printTreeSizes();
    if (this->kain != 0) nNodes += this->kain->printTreeSizes();

    TelePrompter::printSeparator(0, '-');
    println(0," Total number of nodes                   " << setw(18) << nNodes);
    TelePrompter::printSeparator(0, '=', 2);
}

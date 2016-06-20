#include "XCOperator.h"
#include "XCFunctional.h"
#include "FunctionTree.h"
#include "Orbital.h"
#include "Density.h"
#include "Potential.h"
#include "TelePrompter.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

XCOperator::XCOperator(int k,
                       double build_prec,
                       const MultiResolutionAnalysis<3> &mra,
                       XCFunctional &func,
                       OrbitalVector &phi)
        : order(k),
          functional(&func),
          orbitals_0(&phi),
          density_0(func.isSpinSeparated()),
          project(mra),
          grid(mra),
          energy(0.0),
          xcInput(0),
          xcOutput(0) {
    this->potential[0] = new Potential(mra);
    this->potential[1] = new Potential(mra);
    this->potential[2] = new Potential(mra);
}

XCOperator::~XCOperator() {
    if (this->xcInput != 0) MSG_ERROR("XC input not deallocated");
    if (this->xcOutput != 0) MSG_ERROR("XC output not deallocated");
    this->functional = 0;
    this->orbitals_0 = 0;
//    this->gradient_0 = deletePtrArray<Density>(3, &this->gradient_0);
    for (int i = 0; i < 3; i++) {
        if (this->potential[i] == 0) MSG_ERROR("Invalid potential");
        delete this->potential[i];
        this->potential[i] = 0;
    }
}

void XCOperator::calcUnperturbedDensity() {
    if (this->orbitals_0 == 0) MSG_ERROR("Orbitals not initialized");
//    if (this->gradient_0 != 0) MSG_ERROR("Gradient not empty");

    Density &rho = this->density_0;
//    Density **dRho = this->density_0;
    OrbitalVector &phi = *this->orbitals_0;

    this->project.setPrecision(this->apply_prec);

    {
        Timer timer;
        timer.restart();
        this->project(rho, phi);
        double t = timer.getWallTime();
        int n = rho.getNNodes();
        TelePrompter::printTree(0, "XC density", n, t);
    }

    if (this->functional->isGGA()) {
        NOT_IMPLEMENTED_ABORT;
//        Timer timer;
//        timer.restart();
//        dRho = calcDensityGradient(rho);
//        double t = timer.getWallTime();
//        int n = sumNodes<Density>(dRho, 3);
//        TelePrompter::printTree(0, "XC density gradient", n, t);
//        printout(1, endl);
    }
}

/** Compute the required input functions for XCFun. In the case of GGA
 * the spin density gradients are computed. Then densities and gradients
 * are sorted in the input function array in the correct order for XCFun.
 * We define \f$ \gamma_{\alpha\beta} = \nabla\rho_\alpha\cdot\nabla\rho_\beta\f$
 * LDA: \f$ (\rho) \f$
 * GGA: \f$ (\rho, \gamma) \f$
 * Spin LDA: \f$ (\rho_\alpha,\rho_\beta) \f$
 * Spin GGA: \f$ (\rho_\alpha,\rho_\beta, \gamma_{\alpha\alpha},
 * \gamma_{\alpha\beta}, \gamma_{\beta\beta}) \f$
 */
void XCOperator::setupXCInput() {
    if (this->xcInput != 0) MSG_ERROR("XC input not empty");
    Timer timer;
    timer.restart();
    println(2, "Preprocessing");

    int nInp = this->functional->getInputLength();
    this->xcInput = allocPtrArray<FunctionTree<3> >(nInp);

    if (not this->functional->isSpinSeparated()) {
        this->xcInput[0] = &this->density_0.getDensity();
        if (this->functional->isGGA()) {
            NOT_IMPLEMENTED_ABORT;
//            this->xcInput[1] = calcDotProduct<Density>(&this->rho_0[2], &this->rho_0[2]);
        }
    } else {
        this->xcInput[0] = &this->density_0.getDensity(Alpha);
        this->xcInput[1] = &this->density_0.getDensity(Beta);
        if (this->functional->isGGA()) {
            NOT_IMPLEMENTED_ABORT;
//            this->xcInput[2] = calcDotProduct<Density>(&this->rho_0[2], &this->rho_0[2]);
//            this->xcInput[3] = calcDotProduct<Density>(&this->rho_0[2], &this->rho_0[5]);
//            this->xcInput[4] = calcDotProduct<Density>(&this->rho_0[5], &this->rho_0[5]);
        }
    }

    for (int i = 0; i < nInp; i++) {
        if (this->xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }
    double t = timer.getWallTime();
    int n = sumNodes<FunctionTree<3> >(this->xcInput, nInp);
    TelePrompter::printTree(1, "Preprocessing density", n, t);
    printout(2, endl);
}

void XCOperator::clearXCInput() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");

    // These belong to density_0
    this->xcInput[0] = 0;
    if (this->functional->isSpinSeparated()) {
        this->xcInput[1] = 0;
    }

    int nFuncs = this->functional->getInputLength();
    for (int i = 0; i < nFuncs; i++) {
        if (this->xcInput[i] != 0) {
            delete this->xcInput[i];
        }
        this->xcInput[i] = 0;
    }
    delete[] this->xcInput;
    this->xcInput = 0;
}

void XCOperator::setupXCOutput() {
    if (this->xcOutput != 0) MSG_ERROR("XC output not empty");
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");
    if (this->xcInput[0] == 0) MSG_ERROR("XC input not initialized");

    int nOut = this->functional->getOutputLength(this->order);
    this->xcOutput = allocPtrArray<FunctionTree<3> >(nOut);

    // Copy grid from input density
    FunctionTree<3> &rho = *this->xcInput[0];
    for (int i = 0; i < nOut; i++) {
        this->xcOutput[i] = this->grid(rho);
    }
}

void XCOperator::clearXCOutput() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");

    int nFuncs = this->functional->getOutputLength(this->order);
    for (int i = 0; i < nFuncs; i++) {
        if (this->xcOutput[i] != 0) {
            delete this->xcOutput[i];
        }
        this->xcOutput[i] = 0;
    }
    delete[] this->xcOutput;
    this->xcOutput = 0;
}

void XCOperator::evaluateXCFunctional() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");
    if (this->xcOutput == 0) MSG_ERROR("XC input not initialized");

    Timer timer;
    timer.restart();
    println(2, "Evaluating");

    int nInp = this->functional->getInputLength();
    int nOut = this->functional->getOutputLength(this->order);

    MatrixXd inpData, outData;
    compressTreeData(nInp, this->xcInput, inpData);
    this->functional->evaluate(this->order, inpData, outData);
    expandTreeData(nOut, this->xcOutput, outData);

    double t = timer.getWallTime();
    int nNodes = sumNodes<FunctionTree<3> >(this->xcOutput, nOut);
    TelePrompter::printTree(1, "Evaluating XC functional", nNodes, t);
    printout(2, endl);
}

void XCOperator::calcEnergy() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    if (this->xcOutput[0] == 0) MSG_ERROR("Invalid XC output");
    Timer timer;
    timer.restart();

    this->energy = this->xcOutput[0]->integrate();
    double time = timer.getWallTime();
    int nNodes = this->xcOutput[0]->getNNodes();
    TelePrompter::printTree(1, "XC energy", nNodes, time);
}

//Potential** XCOperator::calcGradDotPotDensVec(Potential *pot, Density **dens) {
//    NOT_IMPLEMENTED_ABORT;
//    boost::timer timer;
//    int nNodes;
//    double time;

//    Potential **result = allocPtrArray<Potential>(3);
//    for (int d = 0; d < 3; d++) {
//        if (dens[d] == 0) continue;
//        if (pot == 0) MSG_ERROR("Invalid XC output");

//        timer.restart();
//        Potential potDens;
//        potDens.mult(1.0, *pot, 1.0, *dens[d], 0);
//        time = timer.elapsed();
//        nNodes = potDens.getNNodes();
//        TelePrompter::printTree(2, "Multiply", nNodes, time);

//        timer.restart();
//        result[d] = new Potential;
//        this->derivative->setApplyDir(d);
//        this->derivative->apply(*result[d], potDens);
//        time = timer.elapsed();
//        nNodes = result[d]->getNNodes();
//        TelePrompter::printTree(2, "Gradient", nNodes, time);
//    }
//    return result;
//}

//Potential* calcPotDensVecDotDensVec(Potential *pot, Density **dens_1, Density **dens_2) {
//    NOT_IMPLEMENTED_ABORT;
//}

Orbital* XCOperator::operator() (Orbital &orb) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    switch (orb.getSpin()) {
    case Paired:
        if (this->potential[0] == 0) MSG_ERROR("XC potential not available");
        return (*this->potential[0])(orb);
    case Alpha:
        if (this->potential[1] == 0) MSG_ERROR("Alpha potential not available");
        return (*this->potential[1])(orb);
    case Beta:
        if (this->potential[2] == 0) MSG_ERROR("Beta potential not available");
        return (*this->potential[2])(orb);
    default:
        MSG_ERROR("Invalid spin");
        return 0;
    }
}

Orbital* XCOperator::adjoint(Orbital &orb) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->potential == 0) {
//        setup();
//    }
//    switch (orb.getSpin()) {
//    case Orbital::Paired:
//        if (this->potential[0] == 0) MSG_ERROR("XC potential not available");
//        this->potential[0]->setMaxApplyDepth(this->maxApplyDepth);
//        return this->potential[0]->adjoint(orb);
//    case Orbital::Alpha:
//        if (this->potential[1] == 0) MSG_ERROR("Alpha potential not available");
//        this->potential[1]->setMaxApplyDepth(this->maxApplyDepth);
//        return this->potential[1]->adjoint(orb);
//    case Orbital::Beta:
//        if (this->potential[2] == 0) MSG_ERROR("Beta potential not available");
//        this->potential[2]->setMaxApplyDepth(this->maxApplyDepth);
//        return this->potential[2]->adjoint(orb);
//    default:
//        MSG_ERROR("Invalid spin");
//        return 0;
//    }
}

int XCOperator::printTreeSizes() const {
    NOT_IMPLEMENTED_ABORT;
//    int nNodes = 0;
//    int nTrees = 0;
//    if (this->rho_0 != 0) {
//        int nFuncs = 8;
//        for (int i = 0; i < nFuncs; i++) {
//            if (this->rho_0[i] != 0) {
//                nNodes += this->rho_0[i]->getNNodes();
//                nTrees++;
//            }
//        }
//    }
//    if (this->potential != 0) {
//        if (not this->xcFun->isSpinSeparated()) {
//            if (this->potential[0] != 0) {
//                nNodes += this->potential[0]->getNNodes();
//                nTrees++;
//            }
//        } else {
//            if (this->potential[1] != 0) {
//                nNodes += this->potential[1]->getNNodes();
//                nTrees++;
//            }
//            if (this->potential[2] != 0) {
//                nNodes += this->potential[2]->getNNodes();
//                nTrees++;
//            }
//        }
//    }
//    if (this->xcInput != 0) {
//        int nFuncs = this->xcFun->getInputLength();
//        for (int i = 0; i < nFuncs; i++) {
//            if (this->xcInput[i] != 0) {
//                nNodes += this->xcInput[i]->getNNodes();
//                nTrees++;
//            }
//        }
//    }
//    if (this->xcOutput != 0) {
//        int nFuncs = this->xcFun->getOutputLength();
//        for (int i = 0; i < nFuncs; i++) {
//            if (this->xcOutput[i] != 0) {
//                nNodes += this->xcOutput[i]->getNNodes();
//                nTrees++;
//            }
//        }
//    }
//    println(0, " XCOperator        " << setw(15) << nTrees << setw(25) << nNodes);
//    return nNodes;
}

void XCOperator::compressTreeData(int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");
    if (trees[0] == 0) MSG_ERROR("Invalid input");

    FunctionTree<3> &tree = *trees[0];
    int nCoefs = tree.getTDim()*tree.getKp1_d()*tree.getNEndNodes();
    data = MatrixXd::Zero(nCoefs, nFuncs);

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized input tree");
        VectorXd col_i;
        trees[i]->getEndValues(col_i);
        data.col(i) = col_i;
    }
}

void XCOperator::expandTreeData(int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        trees[i]->setEndValues(col_i);
    }
}

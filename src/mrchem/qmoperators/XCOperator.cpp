#include "XCOperator.h"
#include "XCFunctional.h"
#include "FunctionTree.h"
#include "FunctionNode.h"
#include "Orbital.h"
#include "Density.h"
#include "Potential.h"
#include "TelePrompter.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

XCOperator::XCOperator(int k,
                       double build_prec,
                       XCFunctional &func,
                       OrbitalVector &phi)
        : order(k),
          functional(&func),
          add(-1.0),
          mult(-1.0),
          project(-1.0),
          derivative(-1, *MRA, 0.0, 0.0),
          density_0(func.isSpinSeparated()),
          gradient_0(0),
          orbitals_0(&phi),
          energy(0.0),
          xcInput(0),
          xcOutput(0) {
    this->potential[0] = new Potential();
    this->potential[1] = new Potential();
    this->potential[2] = new Potential();
}

XCOperator::~XCOperator() {
    if (this->xcInput != 0) MSG_ERROR("XC input not deallocated");
    if (this->xcOutput != 0) MSG_ERROR("XC output not deallocated");
    this->functional = 0;
    this->orbitals_0 = 0;
    this->gradient_0 = deletePtrArray<Density>(3, &this->gradient_0);
    for (int i = 0; i < 3; i++) {
        if (this->potential[i] == 0) MSG_ERROR("Invalid potential");
        delete this->potential[i];
        this->potential[i] = 0;
    }
}

void XCOperator::setup(double prec) {
    QMOperator::setup(prec);
    this->add.setPrecision(-1.0);
    this->mult.setPrecision(-1.0);
    this->project.setPrecision(prec);
    this->derivative.setPrecision(prec);
}

void XCOperator::clear() {
    this->add.setPrecision(-1.0);
    this->mult.setPrecision(-1.0);
    this->project.setPrecision(-1.0);
    this->derivative.setPrecision(-1.0);
    QMOperator::clear();
}

void XCOperator::calcDensity() {
    if (this->orbitals_0 == 0) MSG_ERROR("Orbitals not initialized");
    if (this->gradient_0 != 0) MSG_ERROR("Gradient not empty");

    Density &rho = this->density_0;
    OrbitalVector &phi = *this->orbitals_0;

    {
        Timer timer;
        this->project(rho, phi);
        timer.stop();
        double t = timer.getWallTime();
        int n = rho.getNNodes();
        TelePrompter::printTree(0, "XC density", n, t);
    }

    if (this->functional->isGGA()) {
        Timer timer;
        this->gradient_0 = calcDensityGradient(rho);
        timer.stop();
        double t = timer.getWallTime();
        int n = sumNodes<Density>(this->gradient_0, 3);
        TelePrompter::printTree(0, "XC density gradient", n, t);
        printout(1, endl);
    }
}

Density** XCOperator::calcDensityGradient(Density &rho) {
    Density **out = allocPtrArray<Density>(3);

    out[0] = new Density(rho);
    out[1] = new Density(rho);
    out[2] = new Density(rho);

    if (rho.isSpinDensity()) {
        FunctionTree<3> &rho_a = rho.getDensity(Alpha);
        FunctionTreeVector<3> grad_a;
        this->derivative.grad(grad_a, rho_a);
        out[0]->setDensity(Alpha, grad_a[0]);
        out[1]->setDensity(Alpha, grad_a[1]);
        out[2]->setDensity(Alpha, grad_a[2]);

        FunctionTree<3> &rho_b = rho.getDensity(Beta);
        FunctionTreeVector<3> grad_b;
        this->derivative.grad(grad_b, rho_b);
        out[0]->setDensity(Beta, grad_b[0]);
        out[1]->setDensity(Beta, grad_b[1]);
        out[2]->setDensity(Beta, grad_b[2]);
    } else {
        FunctionTree<3> &rho_t = rho.getDensity(Paired);
        FunctionTreeVector<3> grad_t;
        this->derivative.grad(grad_t, rho_t);
        out[0]->setDensity(Paired, grad_t[0]);
        out[1]->setDensity(Paired, grad_t[1]);
        out[2]->setDensity(Paired, grad_t[2]);
    }
    return out;
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
    println(2, "Preprocessing");

    int nInp = this->functional->getInputLength();
    this->xcInput = allocPtrArray<FunctionTree<3> >(nInp);

    if (not this->functional->isSpinSeparated()) {
        this->xcInput[0] = &this->density_0.getDensity(Paired);
        if (this->functional->isGGA()) {
            FunctionTreeVector<3> vec;
            vec.push_back(&this->gradient_0[0]->getDensity(Paired));
            vec.push_back(&this->gradient_0[1]->getDensity(Paired));
            vec.push_back(&this->gradient_0[2]->getDensity(Paired));
            this->xcInput[1] = calcDotProduct(vec, vec);
            vec.clear();
        }
    } else {
        this->xcInput[0] = &this->density_0.getDensity(Alpha);
        this->xcInput[1] = &this->density_0.getDensity(Beta);
        if (this->functional->isGGA()) {
            FunctionTreeVector<3> vec;
            vec.push_back(&this->gradient_0[0]->getDensity(Alpha));
            vec.push_back(&this->gradient_0[1]->getDensity(Alpha));
            vec.push_back(&this->gradient_0[2]->getDensity(Alpha));
            this->xcInput[1] = calcDotProduct(vec, vec);
            vec.clear();

            vec.push_back(&this->gradient_0[0]->getDensity(Beta));
            vec.push_back(&this->gradient_0[1]->getDensity(Beta));
            vec.push_back(&this->gradient_0[2]->getDensity(Beta));
            this->xcInput[2] = calcDotProduct(vec, vec);
            vec.clear();
        }
    }

    for (int i = 0; i < nInp; i++) {
        if (this->xcInput[i] == 0) MSG_ERROR("Invalid XC input");
    }
    timer.stop();
    double t = timer.getWallTime();
    int n = sumNodes<FunctionTree<3> >(this->xcInput, nInp);
    TelePrompter::printTree(0, "XC preprocess xcfun", n, t);
    printout(2, endl);
}

void XCOperator::clearXCInput() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");

    // These belong to density_0
    this->xcInput[0] = 0;
    if (this->functional->isSpinSeparated()) {
        this->xcInput[1] = 0;
    }

    int nInp = this->functional->getInputLength();
    this->xcInput = deletePtrArray<FunctionTree<3> >(nInp, &this->xcInput);
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
        this->xcOutput[i] = new FunctionTree<3>(*MRA);
        this->grid(*this->xcOutput[i], rho);
    }
}

void XCOperator::clearXCOutput() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");

    int nOut = this->functional->getOutputLength(this->order);
    this->xcOutput = deletePtrArray<FunctionTree<3> >(nOut, &this->xcOutput);
}

void XCOperator::evaluateXCFunctional() {
    if (this->xcInput == 0) MSG_ERROR("XC input not initialized");
    if (this->xcOutput == 0) MSG_ERROR("XC input not initialized");

    Timer timer;
    println(2, "Evaluating");

    int nInp = this->functional->getInputLength();
    int nOut = this->functional->getOutputLength(this->order);

#pragma omp parallel firstprivate(nInp, nOut)
{
    	int nNodes = this->xcInput[0]->getNEndNodes();
#pragma omp for schedule(guided)
    	for (int n = 0; n < nNodes; n++) {
            MatrixXd inpData, outData;
            compressNodeData(n, nInp, this->xcInput, inpData);
            this->functional->evaluate(this->order, inpData, outData);
            expandNodeData(n, nOut, this->xcOutput, outData);
        }
}
    for (int i = 0; i < nOut; i++) {
        this->xcOutput[i]->mwTransform(BottomUp);
        this->xcOutput[i]->calcSquareNorm();
    }

    timer.stop();
    double t = timer.getWallTime();
    int n = sumNodes<FunctionTree<3> >(this->xcOutput, nOut);
    TelePrompter::printTree(0, "XC evaluate xcfun", n, t);
    printout(2, endl);
}

void XCOperator::calcEnergy() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    if (this->xcOutput[0] == 0) MSG_ERROR("Invalid XC output");
    Timer timer;
    this->energy = this->xcOutput[0]->integrate();
    timer.stop();
    double time = timer.getWallTime();
    int nNodes = this->xcOutput[0]->getNNodes();
    TelePrompter::printTree(0, "XC energy", nNodes, time);
}

FunctionTree<3>* XCOperator::calcGradDotPotDensVec(FunctionTree<3> &pot,
                                                   FunctionTreeVector<3> &dens) {
    FunctionTreeVector<3> vec;
    for (int d = 0; d < 3; d++) {
        if (dens[d] == 0) MSG_ERROR("Invalid density");

        Timer timer;
        FunctionTree<3> *potDens = new FunctionTree<3>(*MRA);
        this->grid(*potDens, *dens[d]);
        this->mult(*potDens, 1.0, pot, *dens[d], 0);
        vec.push_back(potDens);

        timer.stop();
        double t = timer.getWallTime();
        int n = potDens->getNNodes();
        TelePrompter::printTree(2, "Multiply", n, t);
    }

    Timer timer;
    FunctionTree<3> *result = new FunctionTree<3>(*MRA);
    this->derivative.div(*result, vec);
    vec.clear(true);

    timer.stop();
    double t = timer.getWallTime();
    int n = result->getNNodes();
    TelePrompter::printTree(2, "Gradient", n, t);
    return result;
}

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
    int nNodes = 0;
    nNodes += this->density_0.printTreeSizes();
    this->gradient_0[0]->printTreeSizes();
    this->gradient_0[1]->printTreeSizes();
    this->gradient_0[2]->printTreeSizes();

    if (this->potential != 0) {
        nNodes += this->potential[0]->printTreeSizes();
        nNodes += this->potential[1]->printTreeSizes();
        nNodes += this->potential[2]->printTreeSizes();
    }
    int nInput = 0;
    int inTrees = 0;
    if (this->xcInput != 0) {
        int nFuncs = this->functional->getInputLength();
        for (int i = 0; i < nFuncs; i++) {
            if (this->xcInput[i] != 0) {
                nInput += this->xcInput[i]->getNNodes();
                inTrees++;
            }
        }
    }
    int nOutput = 0;
    int outTrees = 0;
    if (this->xcOutput != 0) {
        int nFuncs = this->functional->getOutputLength(this->order);
        for (int i = 0; i < nFuncs; i++) {
            if (this->xcOutput[i] != 0) {
                nOutput += this->xcOutput[i]->getNNodes();
                outTrees++;
            }
        }
    }
    println(0, " XC input          " << setw(15) << inTrees << setw(25) << nInput);
    println(0, " XC output         " << setw(15) << outTrees << setw(25) << nOutput);

    nNodes += (nInput + nOutput);
    return nNodes;
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

void XCOperator::compressNodeData(int n, int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");
    if (trees[0] == 0) MSG_ERROR("Invalid input");

    FunctionTree<3> &tree = *trees[0];
    int nCoefs = tree.getTDim()*tree.getKp1_d();
    data = MatrixXd::Zero(nCoefs, nFuncs);

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized input tree");
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        VectorXd col_i;
        node.getValues(col_i);
        data.col(i) = col_i;
    }
}

void XCOperator::expandNodeData(int n, int nFuncs, FunctionTree<3> **trees, MatrixXd &data) {
    if (trees == 0) MSG_ERROR("Invalid input");

    for (int i = 0; i < nFuncs; i++) {
        if (trees[i] == 0) MSG_ERROR("Uninitialized output tree " << i);
        VectorXd col_i = data.col(i);
        FunctionNode<3> &node = trees[i]->getEndFuncNode(n);
        node.setValues(col_i);
    }
}

FunctionTree<3>* XCOperator::calcDotProduct(FunctionTreeVector<3> &vec_a,
                                            FunctionTreeVector<3> &vec_b) {
    if (vec_a.size() != vec_b.size()) MSG_ERROR("Invalid input");

    FunctionTreeVector<3> out_vec;
    for (int d = 0; d < vec_a.size(); d++) {
        FunctionTree<3> &tree_a = vec_a.getFunc(d);
        FunctionTree<3> &tree_b = vec_b.getFunc(d);
        FunctionTree<3> *out_d = new FunctionTree<3>(*MRA);
        this->grid(*out_d, tree_a);
        this->grid(*out_d, tree_b);
        this->mult(*out_d, 1.0, tree_a, tree_b, 0);
        out_vec.push_back(out_d);
    }
    FunctionTree<3> *out = new FunctionTree<3>(*MRA);
    this->grid(*out, out_vec);
    this->add(*out, out_vec, 0);

    out_vec.clear(true);
    return out;
}


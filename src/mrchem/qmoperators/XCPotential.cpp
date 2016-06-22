#include "XCPotential.h"
#include "XCFunctional.h"
#include "Potential.h"
#include "TelePrompter.h"

#include "constants.h"

using namespace std;
using namespace Eigen;

XCPotential::XCPotential(double build_prec,
                         const MultiResolutionAnalysis<3> &mra,
                         XCFunctional &func,
                         OrbitalVector &phi)
        : XCOperator(1, build_prec, mra, func, phi) {
}

XCPotential::~XCPotential() {
}

void XCPotential::setup(double prec) {
    this->apply_prec = prec;
    calcUnperturbedDensity();
    setupXCInput();
    setupXCOutput();
    evaluateXCFunctional();
    calcEnergy();
    calcPotential();
    clearXCInput();
    clearXCOutput();
}

void XCPotential::clear() {
    this->apply_prec = -1.0;
    this->energy = 0.0;
    this->density_0.clear();
    this->potential[0]->clear();
    this->potential[1]->clear();
    this->potential[2]->clear();
//    clearPtrArray<Density>(3, this->gradient_0);
}

void XCPotential::calcPotential() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    Timer timer;
    timer.restart();

    if (this->functional->isLDA()) {
        if (not this->functional->isSpinSeparated()) {
            calcPotentialLDA(Paired);
        } else {
            calcPotentialLDA(Alpha);
            calcPotentialLDA(Beta);
        }
    } else if (this->functional->isGGA()) {
        NOT_IMPLEMENTED_ABORT;
    } else {
        MSG_FATAL("Invalid functional type");
    }
    double t = timer.getWallTime();
    int n = -1;
    TelePrompter::printTree(0, "XC potential", n, t);
}

void XCPotential::calcPotentialLDA(int spin) {
    if (spin == Paired) {
        if (this->potential[0] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[0]->real = this->xcOutput[1];
        this->potential[0]->setup(this->apply_prec);
        this->xcOutput[1] = 0;
    } else if (spin == Alpha) {
        if (this->potential[1] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[1]->real = this->xcOutput[1];
        this->potential[1]->setup(this->apply_prec);
        this->xcOutput[1] = 0;
    } else if (spin == Beta) {
        if (this->potential[2] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        this->potential[2]->real = this->xcOutput[2];
        this->potential[2]->setup(this->apply_prec);
        this->xcOutput[2] = 0;
    } else {
        MSG_FATAL("Invalid spin");
    }
}

void XCPotential::calcPotentialGGA(int spin) {
    NOT_IMPLEMENTED_ABORT;
//    FunctionTreeVector<3> xc_funcs;
//    FunctionTreeVector<3> dRho_a;
//    FunctionTreeVector<3> dRho_b;

//    FunctionTree<3> *pot = 0;
//    if (spin == Paired) {
//        xc_funcs.push_back(this->xcOutput[1]);
//        xc_funcs.push_back(this->xcOutput[2]);
//        xc_funcs.push_back(0);
//        dRho_a.push_back(this->gradient_0[0]->total);
//        dRho_a.push_back(this->gradient_0[1]->total);
//        dRho_a.push_back(this->gradient_0[2]->total);
//        dRho_b.push_back(0);
//        dRho_b.push_back(0);
//        dRho_b.push_back(0);

//        pot = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

//        xc_funcs.clear();
//        dRho_a.clear();
//        dRho_b.clear();
//    }
//    if (spin == Alpha) {
//        xc_funcs.push_back(this->xcOutput[1]);
//        xc_funcs.push_back(this->xcOutput[3]);
//        xc_funcs.push_back(this->xcOutput[4]);
//        dRho_a.push_back(this->gradient_0[0]->alpha);
//        dRho_a.push_back(this->gradient_0[1]->alpha);
//        dRho_a.push_back(this->gradient_0[2]->alpha);
//        dRho_b.push_back(this->gradient_0[0]->beta);
//        dRho_b.push_back(this->gradient_0[1]->beta);
//        dRho_b.push_back(this->gradient_0[2]->beta);

//        pot = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

//        xc_funcs.clear();
//        dRho_a.clear();
//        dRho_b.clear();
//    }
//    if (spin == Beta) {
//        xc_funcs.push_back(this->xcOutput[2]);
//        xc_funcs.push_back(this->xcOutput[5]);
//        xc_funcs.push_back(this->xcOutput[4]);
//        dRho_a.push_back(this->gradient_0[0]->beta);
//        dRho_a.push_back(this->gradient_0[1]->beta);
//        dRho_a.push_back(this->gradient_0[2]->beta);
//        dRho_b.push_back(this->gradient_0[0]->alpha);
//        dRho_b.push_back(this->gradient_0[1]->alpha);
//        dRho_b.push_back(this->gradient_0[2]->alpha);

//        pot = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

//        xc_funcs.clear();
//        dRho_a.clear();
//        dRho_b.clear();
//    }
//    return pot;
}

FunctionTree<3>* XCPotential::calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                               FunctionTreeVector<3> &dRho_a,
                                               FunctionTreeVector<3> &dRho_b) {
    NOT_IMPLEMENTED_ABORT;
//    if (xc_funcs[0] == 0) MSG_ERROR("Invalid XC output");

//    FunctionTreeVector<3> funcs;
//    funcs.push_back(1.0, xc_funcs[0]);

//    FunctionTreeVector<3> tmp_1;
//    if (xc_funcs[1] != 0) {
//        tmp_1 = calcGradDotPotDensVec(xc_funcs[1], dRho_a);
//        funcs.push_back(-2.0, tmp_1[0]);
//        funcs.push_back(-2.0, tmp_1[1]);
//        funcs.push_back(-2.0, tmp_1[2]);
//    }

//    FunctionTreeVector<3> tmp_2;
//    if (xc_funcs[2] != 0) {
//        tmp_2 = calcGradDotPotDensVec(xc_funcs[2], dRho_b);
//        funcs.push_back(-1.0, tmp_2[0]);
//        funcs.push_back(-1.0, tmp_2[1]);
//        funcs.push_back(-1.0, tmp_2[2]);
//    }

//    FunctionTree<3> *pot = this->add(funcs);
//    tmp_1.clear(true);
//    tmp_2.clear(true);
//    funcs.clear(false);

//    return pot;
}

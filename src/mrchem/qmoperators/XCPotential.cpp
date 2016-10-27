#include "XCPotential.h"
#include "XCFunctional.h"
#include "Potential.h"
#include "TelePrompter.h"

#include "constants.h"

extern MultiResolutionAnalysis<3> *MRA;

using namespace std;
using namespace Eigen;

void XCPotential::setup(double prec) {
    XCOperator::setup(prec);
    calcDensity();
    setupXCInput();
    setupXCOutput();
    evaluateXCFunctional();
    calcEnergy();
    calcPotential();
    clearXCInput();
    clearXCOutput();
    if (this->potential[0] != 0) this->potential[0]->setup(prec);
    if (this->potential[1] != 0) this->potential[1]->setup(prec);
    if (this->potential[2] != 0) this->potential[2]->setup(prec);
}

void XCPotential::clear() {
    this->energy = 0.0;
    this->density_0.clear();
    if (this->potential[0] != 0) this->potential[0]->clear();
    if (this->potential[1] != 0) this->potential[1]->clear();
    if (this->potential[2] != 0) this->potential[2]->clear();
    if (this->gradient_0 != 0) {
        if (this->gradient_0[0] != 0) this->gradient_0[0]->clear();
        if (this->gradient_0[1] != 0) this->gradient_0[1]->clear();
        if (this->gradient_0[2] != 0) this->gradient_0[2]->clear();
    }
    this->gradient_0 = deletePtrArray<Density>(3, &this->gradient_0);
    XCOperator::clear();
}

void XCPotential::calcPotential() {
    if (this->xcOutput == 0) MSG_ERROR("XC output not initialized");
    Timer timer;
    if (this->functional->isLDA()) {
        if (not this->functional->isSpinSeparated()) {
            calcPotentialLDA(Paired);
        } else {
            calcPotentialLDA(Alpha);
            calcPotentialLDA(Beta);
        }
    } else if (this->functional->isGGA()) {
        if (not this->functional->isSpinSeparated()) {
            calcPotentialGGA(Paired);
        } else {
            calcPotentialGGA(Alpha);
            calcPotentialGGA(Beta);
        }
    } else {
        MSG_FATAL("Invalid functional type");
    }
    timer.stop();
    double t = timer.getWallTime();
    int n = 0;
    if (this->potential[0] != 0) n += this->potential[0]->getNNodes();
    if (this->potential[1] != 0) n += this->potential[1]->getNNodes();
    if (this->potential[2] != 0) n += this->potential[2]->getNNodes();
    TelePrompter::printTree(0, "XC potential", n, t);
}

void XCPotential::calcPotentialLDA(int spin) {
    if (spin == Paired) {
        if (this->potential[0] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[0]->real = this->xcOutput[1];
        this->xcOutput[1] = 0;
    } else if (spin == Alpha) {
        if (this->potential[1] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        this->potential[1]->real = this->xcOutput[1];
        this->xcOutput[1] = 0;
    } else if (spin == Beta) {
        if (this->potential[2] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        this->potential[2]->real = this->xcOutput[2];
        this->xcOutput[2] = 0;
    } else {
        MSG_FATAL("Invalid spin");
    }
}

void XCPotential::calcPotentialGGA(int spin) {
    FunctionTreeVector<3> xc_funcs;
    FunctionTreeVector<3> dRho_a;
    FunctionTreeVector<3> dRho_b;

    if (spin == Paired) {
        if (this->potential[0] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[1]);
        xc_funcs.push_back(this->xcOutput[2]);
        xc_funcs.push_back(0);
        dRho_a.push_back(&this->gradient_0[0]->getDensity(Paired));
        dRho_a.push_back(&this->gradient_0[1]->getDensity(Paired));
        dRho_a.push_back(&this->gradient_0[2]->getDensity(Paired));
        dRho_b.push_back(0);
        dRho_b.push_back(0);
        dRho_b.push_back(0);

        this->potential[0]->real = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
    if (spin == Alpha) {
        if (this->potential[1] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[1] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[3] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[4] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[1]);
        xc_funcs.push_back(this->xcOutput[3]);
        xc_funcs.push_back(this->xcOutput[4]);
        dRho_a.push_back(&this->gradient_0[0]->getDensity(Alpha));
        dRho_a.push_back(&this->gradient_0[1]->getDensity(Alpha));
        dRho_a.push_back(&this->gradient_0[2]->getDensity(Alpha));
        dRho_b.push_back(&this->gradient_0[0]->getDensity(Beta));
        dRho_b.push_back(&this->gradient_0[1]->getDensity(Beta));
        dRho_b.push_back(&this->gradient_0[2]->getDensity(Beta));

        this->potential[1]->real = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
    if (spin == Beta) {
        if (this->potential[2] == 0) MSG_ERROR("Invalid XC potential");
        if (this->xcOutput[2] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[4] == 0) MSG_ERROR("Invalid XC output");
        if (this->xcOutput[5] == 0) MSG_ERROR("Invalid XC output");
        xc_funcs.push_back(this->xcOutput[2]);
        xc_funcs.push_back(this->xcOutput[5]);
        xc_funcs.push_back(this->xcOutput[4]);
        dRho_a.push_back(&this->gradient_0[0]->getDensity(Beta));
        dRho_a.push_back(&this->gradient_0[1]->getDensity(Beta));
        dRho_a.push_back(&this->gradient_0[2]->getDensity(Beta));
        dRho_b.push_back(&this->gradient_0[0]->getDensity(Alpha));
        dRho_b.push_back(&this->gradient_0[1]->getDensity(Alpha));
        dRho_b.push_back(&this->gradient_0[2]->getDensity(Alpha));

        this->potential[2]->real = calcPotentialGGA(xc_funcs, dRho_a, dRho_b);

        xc_funcs.clear();
        dRho_a.clear();
        dRho_b.clear();
    }
}

FunctionTree<3>* XCPotential::calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                               FunctionTreeVector<3> &dRho_a,
                                               FunctionTreeVector<3> &dRho_b) {
    if (xc_funcs[0] == 0) MSG_ERROR("Invalid XC output");

    FunctionTreeVector<3> funcs;
    funcs.push_back(1.0, xc_funcs[0]);

    FunctionTree<3> *tmp_1 = 0;
    if (xc_funcs[1] != 0) {
        tmp_1 = calcGradDotPotDensVec(*xc_funcs[1], dRho_a);
        funcs.push_back(-2.0, tmp_1);
    }

    FunctionTree<3> *tmp_2 = 0;
    if (xc_funcs[2] != 0) {
        tmp_2 = calcGradDotPotDensVec(*xc_funcs[2], dRho_b);
        funcs.push_back(-1.0, tmp_2);
    }

    FunctionTree<3> *pot = new FunctionTree<3>(*MRA);
    this->grid(*pot, funcs);
    this->add(*pot, funcs, 0);
    funcs.clear(false);

    if (tmp_1 != 0) delete tmp_1;
    if (tmp_2 != 0) delete tmp_2;

    return pot;
}

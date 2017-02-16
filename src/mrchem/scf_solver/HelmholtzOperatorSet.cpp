#include "HelmholtzOperatorSet.h"
#include "OrbitalVector.h"
#include "OperatorTree.h"
#include "Timer.h"
#include "eigen_disable_warnings.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

HelmholtzOperatorSet::HelmholtzOperatorSet(double build, double thrs)
    : threshold(thrs),
      build_prec(build),
      apply(-1.0, MRA->getMaxScale(), true) {
}

void HelmholtzOperatorSet::initialize(const VectorXd &energies) {
    TelePrompter::printHeader(0, "Initializing Helmholtz Operators");

    Timer timer;
    this->operIdx.clear();
    this->lambda.clear();
    for (int i = 0; i < energies.size(); i++) {
        double energy = energies(i);
 
        if (energy > 0.0) {
            energy = -0.5;
        }
        int idx = initHelmholtzOperator(energy);
        this->operIdx.push_back(idx);
        double mu = getOperator(i).getMu();
        this->lambda.push_back(-0.5*mu*mu);
    }
    clearUnused();
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}

int HelmholtzOperatorSet::initHelmholtzOperator(double energy) {
    if (energy > 0.0) MSG_ERROR("Complex Helmholtz not available");
    double mu = sqrt(-2.0*energy);
    for (int j = 0; j < this->operators.size(); j++) {
        double mu_j = this->operators[j]->getMu();
        double muDiff = mu - mu_j;
        double relDiff = muDiff/mu;
        if (fabs(relDiff) < this->threshold) {
            double l = -0.5*mu_j*mu_j;
            TelePrompter::printDouble(0, "Re-using operator with lambda", l);
            return j;
        }
    }
    TelePrompter::printDouble(0, "Creating operator with lambda", energy);

    HelmholtzOperator *oper = new HelmholtzOperator(*MRA, mu, this->build_prec);
    this->operators.push_back(oper);

    return this->operators.size() - 1;
}

void HelmholtzOperatorSet::clear() {
    int nOper = this->operators.size();
    for (int j = 0; j < nOper; j++) {
        if (this->operators[j] != 0) {
            delete this->operators[j];
            this->operators[j] = 0;
        }
    }
    this->operators.clear();
    this->operIdx.clear();
    this->lambda.clear();
}

void HelmholtzOperatorSet::clearUnused() {
    vector<HelmholtzOperator *> tmp;
    int nIdx = this->operIdx.size();
    int nOper = this->operators.size();
    VectorXd nUsed = VectorXd::Zero(nOper);
    for (int i = 0; i < nIdx; i++) {
        int n = this->operIdx[i];
        nUsed(n) = nUsed(n) + 1;
    }
    for (int i = 0; i < nOper; i++) {
        if (nUsed(i) == 0) {
            delete this->operators[i];
        } else {
            tmp.push_back(this->operators[i]);
            int newIdx = tmp.size() - 1;
            for (int j = 0; j < nIdx; j++) {
                if (operIdx[j] == i) operIdx[j] = newIdx;
            }
        }
        this->operators[i] = 0;
    }
    this->operators.clear();
    for (int i = 0; i < tmp.size(); i++) {
        this->operators.push_back(tmp[i]);
        tmp[i] = 0;
    }
}

HelmholtzOperator& HelmholtzOperatorSet::getOperator(int i) {
    int idx = this->operIdx[i];
    if (idx < 0 or idx >= operators.size()) MSG_ERROR("Invalid operator index");
    HelmholtzOperator *oper = this->operators[idx];
    if (oper == 0) MSG_ERROR("Operator null pointer");
    return *oper;
}

VectorXd HelmholtzOperatorSet::getLambda() const {
    int nLambda = this->lambda.size();
    VectorXd L = VectorXd::Zero(nLambda);
    for (int i = 0; i < nLambda; i++) {
        L(i) = this->lambda[i];
    }
    return L;
}

/** Prints the number of trees and nodes kept in the operator set */
/*
int HelmholtzOperatorSet::printTreeSizes() const {
    int totNodes = 0;
    int totTrees = 0;
    int nOperators = this->operators.size();
    for (int i = 0; i < nOperators; i++) {
        int nTrees = this->operators[i]->size();
        for (int j = 0; j < nTrees; j++) {
            totNodes += this->operators[i]->getComponent(j).getNNodes();
        }
        totTrees += nTrees;
    }
    println(0, " HelmholtzOperators" << setw(15) << totTrees << setw(25) << totNodes);
    return totNodes;
}
*/

void HelmholtzOperatorSet::operator()(int i, Orbital &out, Orbital &inp) {
    if (out.hasReal()) MSG_ERROR("Orbital not empty");
    if (out.hasImag()) MSG_ERROR("Orbital not empty");

    HelmholtzOperator &H_i = getOperator(i);

    if (inp.hasReal()) {
        out.allocReal();
        this->apply(out.real(), H_i, inp.real());
    }
    if (inp.hasImag()) {
        out.allocImag();
        this->apply(out.imag(), H_i, inp.imag());
    }
}

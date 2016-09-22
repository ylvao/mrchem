#include "FockOperator.h"
#include "KineticOperator.h"
#include "NuclearPotential.h"
#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "XCOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

FockOperator::FockOperator(const MultiResolutionAnalysis<3> &mra,
                           KineticOperator *t,
                           NuclearPotential *v,
                           CoulombOperator *j,
                           ExchangeOperator *k,
                           XCOperator *xc)
        : QMOperator(mra),
          add(mra, -1.0),
          T(t),
          V(v),
          J(j),
          K(k),
          XC(xc) {
}

FockOperator::~FockOperator() {
    this->T = 0;
    this->V = 0;
    this->J = 0;
    this->K = 0;
    this->XC = 0;
}

void FockOperator::setup(double prec) {
    Timer timer;
    TelePrompter::printHeader(0, "Setting up Fock operator");
    QMOperator::setup(prec);
    this->add.setPrecision(prec);
    if (this->T != 0) this->T->setup(prec);
    if (this->V != 0) this->V->setup(prec);
    if (this->J != 0) this->J->setup(prec);
    if (this->K != 0) this->K->setup(prec);
    if (this->XC != 0) this->XC->setup(prec);
    for (int i = 0; i < getNPerturbations(); i++) {
        getPerturbationOperator(i).setup(prec);
    }
    timer.stop();
    TelePrompter::printFooter(0, timer, 2);
}

void FockOperator::clear() {
    this->add.setPrecision(-1.0);
    if (this->T != 0) this->T->clear();
    if (this->V != 0) this->V->clear();
    if (this->J != 0) this->J->clear();
    if (this->K != 0) this->K->clear();
    if (this->XC != 0) this->XC->clear();
    for (int i = 0; i < getNPerturbations(); i++) {
        getPerturbationOperator(i).clear();
    }
    QMOperator::clear();
}

void FockOperator::rotate(MatrixXd &U) {
    if (this->T != 0) this->T->rotate(U);
    if (this->V != 0) this->V->rotate(U);
    if (this->J != 0) this->J->rotate(U);
    if (this->K != 0) this->K->rotate(U);
    if (this->XC != 0) this->XC->rotate(U);
    for (int i = 0; i < getNPerturbations(); i++) {
        getPerturbationOperator(i).rotate(U);
    }
}

Orbital* FockOperator::operator() (Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
//    Timer timer;
//    FunctionTreeVector orbs;
//    if (this->T != 0) orbs.push_back((*this->T)(orb_p));
//    if (this->V != 0) orbs.push_back((*this->V)(orb_p));
//    if (this->J != 0) orbs.push_back((*this->J)(orb_p));
//    if (this->K != 0) orbs.push_back((*this->K)(orb_p));
//    if (this->XC != 0) orbs.push_back((*this->XC)(orb_p));
//    for (int i = 0; i < getNPerturbations(); i++) {
//        QMOperator &h1 = getPerturbationOperator(i);
//        orbs.push_back(h1(orb_p));
//    }

//    Orbital *result = new Orbital(orb_p);
//    ->add(orbs, 0);
//    double time = timer.elapsed();
//    int nNodes = result->getNNodes();
//    TelePrompter::printTree(1, "Sum Fock operator", nNodes, time);
//    for (int n = 0; n < orbs.size(); n++) {
//        delete orbs[n];
//        orbs[n] = 0;
//    }
//    return result;
}

Orbital* FockOperator::adjoint(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
//    vector<FunctionTree<3> *> orbs;
//    if (this->T != 0) orbs.push_back((*this->T).adjoint(orb_p));
//    if (this->V != 0) orbs.push_back((*this->V).adjoint(orb_p));
//    if (this->J != 0) orbs.push_back((*this->J).adjoint(orb_p));
//    if (this->K != 0) orbs.push_back((*this->K).adjoint(orb_p));
//    if (this->XC != 0) orbs.push_back((*this->XC).adjoint(orb_p));
//    for (int i = 0; i < getNPerturbations(); i++) {
//        QMOperator &h1 = getPerturbationOperator(i);
//        orbs.push_back(h1.adjoint(orb_p));
//    }

//    Orbital *result = new Orbital(orb_p);
//    result->add(orbs, 0);
//    double time = timer.elapsed();
//    int nNodes = result->getNNodes();
//    TelePrompter::printTree(1, "Sum Fock operator", nNodes, time);
//    for (int n = 0; n < orbs.size(); n++) {
//        delete orbs[n];
//        orbs[n] = 0;
//    }
//    return result;
}

double FockOperator::operator() (Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->T != 0) result += (*this->T)(orb_i, orb_j);
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h1 = getPerturbationOperator(i);
        result += h1(orb_i, orb_j);
    }
    return result;
}

double FockOperator::adjoint(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->T != 0) result += (*this->T).adjoint(orb_i, orb_j);
    if (this->V != 0) result += (*this->V).adjoint(orb_i, orb_j);
    if (this->J != 0) result += (*this->J).adjoint(orb_i, orb_j);
    if (this->K != 0) result += (*this->K).adjoint(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC).adjoint(orb_i, orb_j);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h1 = getPerturbationOperator(i);
        result += h1.adjoint(orb_i, orb_j);
    }
    return result;
}

MatrixXd FockOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);
    if (this->T != 0) result += (*this->T)(i_orbs, j_orbs);
    if (this->V != 0) result += (*this->V)(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J)(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K)(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC)(i_orbs, j_orbs);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h1 = getPerturbationOperator(i);
        result += h1(i_orbs, j_orbs);
    }
    return result;
}

MatrixXd FockOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXd result = MatrixXd::Zero(Ni,Nj);
    if (this->T != 0) result += (*this->T).adjoint(i_orbs, j_orbs);
    if (this->V != 0) result += (*this->V).adjoint(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J).adjoint(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K).adjoint(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC).adjoint(i_orbs, j_orbs);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h1 = getPerturbationOperator(i);
        result += h1.adjoint(i_orbs, j_orbs);
    }
    return result;
}

Orbital* FockOperator::applyKinetic(Orbital &orb_p) {
    return (*this->T)(orb_p);
}

double FockOperator::applyKinetic(Orbital &orb_i, Orbital &orb_j) {
    return (*this->T)(orb_i,orb_j);
}

MatrixXd FockOperator::applyKinetic(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    return (*this->T)(i_orbs,j_orbs);
}

Orbital* FockOperator::applyAdjointKinetic(Orbital &orb_p) {
    return (*this->T).adjoint(orb_p);
}

double FockOperator::applyAdjointKinetic(Orbital &orb_i, Orbital &orb_j) {
    return (*this->T).adjoint(orb_i,orb_j);
}

MatrixXd FockOperator::applyAdjointKinetic(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    return (*this->T).adjoint(i_orbs,j_orbs);
}

Orbital* FockOperator::applyPotential(Orbital &orb_p) {
    vector<double> coefs;
    vector<Orbital *> orbs;
    if (this->V != 0) {
        coefs.push_back(1.0);
        orbs.push_back((*this->V)(orb_p));
    }
    if (this->J != 0) {
        coefs.push_back(1.0);
        orbs.push_back((*this->J)(orb_p));
    }
    if (this->K != 0) {
        coefs.push_back(1.0);
        orbs.push_back((*this->K)(orb_p));
    }
    if (this->XC != 0) {
        coefs.push_back(1.0);
        orbs.push_back((*this->XC)(orb_p));
    }

    Timer timer;
    Orbital *result = new Orbital(orb_p);
    this->add(*result, coefs, orbs, false);
    timer.stop();
    double time = timer.getWallTime();
    int nNodes = result->getNNodes();
    TelePrompter::printTree(1, "Sum potential operator", nNodes, time);

    for (int n = 0; n < orbs.size(); n++) {
        delete orbs[n];
        orbs[n] = 0;
    }
    return result;
}

double FockOperator::applyPotential(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->V != 0) result += (*this->V)(orb_i, orb_j);
    if (this->J != 0) result += (*this->J)(orb_i, orb_j);
    if (this->K != 0) result += (*this->K)(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC)(orb_i, orb_j);
    return result;
}

MatrixXd FockOperator::applyPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    if (this->V != 0) result += (*this->V)(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J)(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K)(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC)(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyAdjointPotential(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
//    vector<FunctionTree<3> *> orbs;

//    if (this->V != 0) orbs.push_back((*this->V).adjoint(orb_p));
//    if (this->J != 0) orbs.push_back((*this->J).adjoint(orb_p));
//    if (this->K != 0) orbs.push_back((*this->K).adjoint(orb_p));
//    if (this->XC != 0) orbs.push_back((*this->XC).adjoint(orb_p));

//    Orbital *result = new Orbital(orb_p);
//    result->add(orbs, 0);
//    double time = timer.elapsed();
//    int nNodes = result->getNNodes();
//    TelePrompter::printTree(1, "Sum adjoint potential operator", nNodes, time);
//    for (int n = 0; n < orbs.size(); n++) {
//        delete orbs[n];
//        orbs[n] = 0;
//    }
//    return result;
}

double FockOperator::applyAdjointPotential(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    if (this->V != 0) result += (*this->V).adjoint(orb_i, orb_j);
    if (this->J != 0) result += (*this->J).adjoint(orb_i, orb_j);
    if (this->K != 0) result += (*this->K).adjoint(orb_i, orb_j);
    if (this->XC != 0) result += (*this->XC).adjoint(orb_i, orb_j);
    return result;
}

MatrixXd FockOperator::applyAdjointPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    if (this->V != 0) result += (*this->V).adjoint(i_orbs, j_orbs);
    if (this->J != 0) result += (*this->J).adjoint(i_orbs, j_orbs);
    if (this->K != 0) result += (*this->K).adjoint(i_orbs, j_orbs);
    if (this->XC != 0) result += (*this->XC).adjoint(i_orbs, j_orbs);
    return result;
}

Orbital* FockOperator::applyPerturbations(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
//    vector<FunctionTree<3> *> orbs;

//    for (int i = 0; i < getNPerturbations(); i++) {
//        QMOperator &h = getPerturbationOperator(i);
//        orbs.push_back(h(orb_p));
//    }

//    Orbital *result = new Orbital(orb_p);
//    result->add(orbs, 0);
//    double time = timer.elapsed();
//    int nNodes = result->getNNodes();
//    TelePrompter::printTree(1, "Sum perturbation operator", nNodes, time);
//    for (int n = 0; n < orbs.size(); n++) {
//        delete orbs[n];
//        orbs[n] = 0;
//    }
//    return result;
}

double FockOperator::applyPerturbations(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h = getPerturbationOperator(i);
        result += h(orb_i, orb_j);
    }
    return result;
}

MatrixXd FockOperator::applyPerturbations(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h = getPerturbationOperator(i);
        result += h(i_orbs, j_orbs);
    }
    return result;
}

Orbital* FockOperator::applyAdjointPerturbations(Orbital &orb_p) {
    NOT_IMPLEMENTED_ABORT;
//    vector<FunctionTree<3> *> orbs;

//    for (int i = 0; i < getNPerturbations(); i++) {
//        QMOperator &h = getPerturbationOperator(i);
//        orbs.push_back(h.adjoint(orb_p));
//    }

//    Orbital *result = new Orbital(orb_p);
//    result->add(orbs, 0);
//    double time = timer.elapsed();
//    int nNodes = result->getNNodes();
//    TelePrompter::printTree(1, "Sum perturbation operator", nNodes, time);
//    for (int n = 0; n < orbs.size(); n++) {
//        delete orbs[n];
//        orbs[n] = 0;
//    }
//    return result;
}

double FockOperator::applyAdjointPerturbations(Orbital &orb_i, Orbital &orb_j) {
    double result = 0.0;
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h = getPerturbationOperator(i);
        result += h(orb_i, orb_j);
    }
    return result;
}

MatrixXd FockOperator::applyAdjointPerturbations(OrbitalVector &i_orbs, OrbitalVector &j_orbs) {
    int nOrbs = i_orbs.size();
    MatrixXd result = MatrixXd::Zero(nOrbs,nOrbs);
    for (int i = 0; i < getNPerturbations(); i++) {
        QMOperator &h = getPerturbationOperator(i);
        result += h(i_orbs, j_orbs);
    }
    return result;
}

/** Prints the number of trees and nodes kept in the Fock operator */
int FockOperator::printTreeSizes() const {
    int nNodes = 0;
    if (this->V != 0) nNodes += this->V->printTreeSizes();
    if (this->J != 0) nNodes += this->J->printTreeSizes();
    if (this->K != 0) nNodes += this->K->printTreeSizes();
    if (this->XC != 0) nNodes += this->XC->printTreeSizes();

    for (int i = 0; i < getNPerturbations(); i++) {
        const QMOperator &h1 = getPerturbationOperator(i);
        nNodes += h1.printTreeSizes();
    }
    return nNodes;
}

#include "CoulombOperator.h"
#include "FunctionTree.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

CoulombOperator::CoulombOperator(double prec, OrbitalVector &phi)
    : QMOperator(MRA->getMaxScale()),
      poisson(*MRA, prec),
      project(-1.0),
      density(Paired),
      potential(),
      orbitals(&phi) {
}

void CoulombOperator::setup(double prec) {
    QMOperator::setup(prec);
    this->project.setPrecision(this->apply_prec);
}

void CoulombOperator::clear() {
    this->project.setPrecision(-1.0);
    QMOperator::clear();
}

int CoulombOperator::printTreeSizes() const {
    int nNodes = 0;
    nNodes += this->density.printTreeSizes();
    nNodes += this->potential.printTreeSizes();
    return nNodes;
}

Orbital* CoulombOperator::operator() (Orbital &orb) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    return this->potential(orb);
}

Orbital* CoulombOperator::adjoint(Orbital &orb) {
    NOT_IMPLEMENTED_ABORT;
}

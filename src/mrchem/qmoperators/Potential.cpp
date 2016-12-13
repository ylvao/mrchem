#include "Potential.h"
#include "FunctionTreeVector.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;
using namespace Eigen;

Potential::Potential()
    : ComplexFunction<3>(0, 0),
      QMOperator(MRA->getMaxScale()) {
}

void Potential::setup(double prec) {
    if (this->hasReal()) MSG_ERROR("Operator not cleared");
    if (this->hasImag()) MSG_ERROR("Operator not cleared");
    QMOperator::setup(prec);
}

void Potential::clear() {
    clearReal(true);
    clearImag(true);
    QMOperator::clear();
}

Orbital* Potential::operator() (Orbital &phi_p) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    OrbitalMultiplier mult(this->apply_prec, this->max_scale);

    Potential &V = *this;
    Orbital *Vphi_p = new Orbital(phi_p);

    Timer timer;
    mult(*Vphi_p, V, phi_p);
    timer.stop();

    int n = Vphi_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied potential", n, t);

    return Vphi_p;
}

Orbital* Potential::adjoint(Orbital &orb) {
    NOT_IMPLEMENTED_ABORT;
}


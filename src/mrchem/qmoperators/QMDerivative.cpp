#include "QMDerivative.h"
#include "MWDerivative.h"
#include "DerivativeOperator.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

QMDerivative::QMDerivative(int dir, DerivativeOperator<3> &d)
    : QMOperator(MRA->getMaxScale()),
      apply_dir(dir),
      derivative(&d) {
}

Orbital* QMDerivative::operator()(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Timer timer;
    Orbital *dPhi = new Orbital(phi);
    calcRealPart(*dPhi, phi);
    calcImagPart(*dPhi, phi);
    timer.stop();

    int n = dPhi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied QM derivative", n, t);

    return dPhi;
}

Orbital* QMDerivative::adjoint(Orbital &phi) {
    NOT_IMPLEMENTED_ABORT;
    return 0;
}

void QMDerivative::calcRealPart(Orbital &dPhi, Orbital &phi) {
    if (this->derivative == 0) MSG_ERROR("No derivative operator");
    MWDerivative<3> apply(this->max_scale);

    if (this->isReal() and phi.hasReal()) {
        dPhi.allocReal();
        apply(dPhi.real(), *this->derivative, phi.real(), this->apply_dir);
    }
    if (this->isImag() and phi.hasImag()) {
        dPhi.allocReal();
        apply(dPhi.real(), *this->derivative, phi.imag(), this->apply_dir);
        dPhi.real() *= -1.0;
    }
}

void QMDerivative::calcImagPart(Orbital &dPhi, Orbital &phi) {
    if (this->derivative == 0) MSG_ERROR("No derivative operator");
    MWDerivative<3> apply(this->max_scale);

    if (this->isImag() and phi.hasReal()) {
        dPhi.allocImag();
        apply(dPhi.imag(), *this->derivative, phi.real(), this->apply_dir);
    }
    if (this->isReal() and phi.hasImag()) {
        dPhi.allocImag();
        apply(dPhi.imag(), *this->derivative, phi.imag(), this->apply_dir);
    }
}

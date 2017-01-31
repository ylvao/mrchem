#include "QMDerivative.h"
#include "MWDerivative.h"
#include "ABGVOperator.h"
#include "PHOperator.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

QMDerivative::QMDerivative(int dir)
    : QMOperator(MRA->getMaxScale()),
      apply_dir(dir),
      //diff_oper(*MRA, 1) {
      diff_oper(*MRA, 0.0, 0.0) {
}

Orbital* QMDerivative::operator()(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    MWDerivative<3> apply(this->max_scale);

    DerivativeOperator<3> &D = this->diff_oper;
    Orbital *dPhi = new Orbital(phi);

    Timer timer;
    if (phi.hasReal()) {
        dPhi->allocReal();
        apply(dPhi->real(), D, phi.real(), this->apply_dir);
    }
    if (phi.hasImag()) {
        dPhi->allocImag();
        apply(dPhi->imag(), D, phi.imag(), this->apply_dir);
    }
    timer.stop();

    int n = dPhi->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(2, "Applied QM derivative", n, t);

    return dPhi;
}

Orbital* QMDerivative::adjoint(Orbital &phi) {
    Orbital *dPhi = (*this)(phi);
    return dPhi;
}

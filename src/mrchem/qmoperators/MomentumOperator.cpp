#include "MomentumOperator.h"
#include "MWDerivative.h"
#include "ABGVOperator.h"
#include "PHOperator.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

MomentumOperator::MomentumOperator(int dir, double prec)
    : QMOperator(MRA->getMaxScale()),
      apply_dir(dir),
      diff_oper(*MRA, 0.5, 0.5) {
}

void MomentumOperator::setup(double prec) {
    QMOperator::setup(prec);
}

void MomentumOperator::clear() {
    QMOperator::clear();
}

Orbital* MomentumOperator::operator() (Orbital &orb_p) {
    Timer timer;
    MWDerivative<3> apply(this->max_scale);
    Orbital *dOrb_p = new Orbital(orb_p);
    if (orb_p.hasReal()) {
        dOrb_p->allocImag();
        apply(dOrb_p->im(), this->diff_oper, orb_p.re(), this->apply_dir);
    }
    if (orb_p.hasImag()) {
        dOrb_p->allocReal();
        apply(dOrb_p->re(), this->diff_oper, orb_p.im(), this->apply_dir);
        dOrb_p->re() *= -1.0;
    }
    timer.stop();
    int n = dOrb_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(2, "Applied momentum", n, t);
    return dOrb_p;
}

Orbital* MomentumOperator::adjoint(Orbital &orb_p) {
    Orbital *dOrb_p = (*this)(orb_p);
    return dOrb_p;
}

int MomentumOperator::printTreeSizes() const {
    println(0, " MomentumOperator  " << setw(15) << 0 << setw(25) << 0);
    return 0;
}

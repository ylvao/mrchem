#include "MomentumOperator.h"
#include "Orbital.h"
#include "Timer.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

using namespace std;

MomentumOperator::MomentumOperator(int dir, double build_prec)
        : apply_dir(dir),
          derivative(*MRA, 0.0, 0.0),
          apply(-1.0, MRA->getMaxScale()) {
}

void MomentumOperator::setup(double prec) {
    QMOperator::setup(prec);
    this->apply.setPrecision(prec);
}

void MomentumOperator::clear() {
    this->apply.setPrecision(-1.0);
    QMOperator::clear();
}

Orbital* MomentumOperator::operator() (Orbital &orb_p) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    Timer timer;
    Orbital *dOrb_p = new Orbital(orb_p);
    if (orb_p.hasReal()) {
        dOrb_p->allocImag();
        this->grid(dOrb_p->im(), orb_p.re());
        this->apply(dOrb_p->im(), this->derivative, orb_p.re(), 0, this->apply_dir);
    }
    if (orb_p.hasImag()) {
        dOrb_p->allocReal();
        this->grid(dOrb_p->re(), orb_p.im());
        this->apply(dOrb_p->re(), this->derivative, orb_p.im(), 0, this->apply_dir);
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

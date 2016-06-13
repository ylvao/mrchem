#include "MomentumOperator.h"
#include "DerivativeOperator.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;

MomentumOperator::MomentumOperator(int dir, DerivativeOperator<3> &d_oper)
        : QMOperator(true),
          apply_dir(dir),
          D(&d_oper) {
}

/* Derivative operator does not belong to this class and is not deleted */
MomentumOperator::~MomentumOperator() {
    this->D = 0;
}

Orbital* MomentumOperator::operator() (Orbital &orb_p) {
    Timer timer;
    timer.restart();
    Orbital *dOrb_p = new Orbital(orb_p);
    this->D->setApplyDir(this->apply_dir);
    if (orb_p.real != 0) {
        dOrb_p->imag = (*this->D)(*orb_p.real);
    }
    if (orb_p.imag != 0) {
        dOrb_p->real = (*this->D)(*orb_p.imag);
        *dOrb_p->real *= -1.0;
    }
    timer.stop();
    int n = dOrb_p->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Applied momentum", n, t);
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

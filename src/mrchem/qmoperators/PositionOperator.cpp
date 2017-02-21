#include "PositionOperator.h"
#include "MWProjector.h"

void QMPosition::setPosition(double *out, const double *inp) {
    if (inp != 0) {
        out[0] = inp[0];
        out[1] = inp[1];
        out[2] = inp[2];
    } else {
        out[0] = 0.0;
        out[1] = 0.0;
        out[2] = 0.0;
    }
}

void QMPosition::setFunction(int d, const double *o) {
    this->func = [d, o] (const double *r) -> double { return r[d]-o[d]; };
}

void QMPosition::setup(double prec) {
    if (IS_EQUAL(prec, this->apply_prec)) return;

    setApplyPrec(prec);
    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    MWProjector<3> project(this->apply_prec, this->max_scale);

    Timer timer;
    this->allocReal();
    project(this->real(), this->func);
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Position operator", n, t);
}

void QMPosition::clear() {
    clearReal(true);
    clearImag(true);
    clearApplyPrec();
}

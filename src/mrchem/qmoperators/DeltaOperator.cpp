#include "DeltaOperator.h"
#include "GridGenerator.h"
#include "MWProjector.h"

void DeltaOperator::setPosition(double *out, const double *inp) {
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

void DeltaOperator::setFunction(double beta, const double *pos) {
    double alpha = pow(beta/pi, 3.0/2.0);
    this->func.setCoef(alpha);
    this->func.setExp(beta);
    this->func.setPos(pos);
}

void DeltaOperator::setup(double prec) {
    if (IS_EQUAL(prec, this->apply_prec)) return;

    setApplyPrec(prec);
    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    GridGenerator<3> grid(this->max_scale);
    MWProjector<3> project(this->apply_prec, this->max_scale);

    Timer timer;
    this->allocReal();
    grid(this->real(), this->func);
    project(this->real(), this->func);
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    TelePrompter::printTree(1, "Position operator", n, t);
}

void DeltaOperator::clear() {
    clearReal(true);
    clearImag(true);
    clearApplyPrec();
}

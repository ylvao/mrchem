#ifndef DELTAOPERATOR_H
#define DELTAOPERATOR_H

#include "QMOperator.h"

class DeltaOperator : public QMOperator {
public:
    DeltaOperator(const double *o) { }
    virtual ~SpinOperator() { }

    virtual void setup(double prec) { NOT_IMPLEMENTED_ABORT; }
    virtual void clear() { NOT_IMPLEMENTED_ABORT; }

    virtual Orbital* operator()(Orbital &orb_p) { NOT_IMPLEMENTED_ABORT; }
    virtual Orbital* adjoint(Orbital &orb_p) { NOT_IMPLEMENTED_ABORT; }

    using QMOperator::operator();
    using QMOperator::adjoint;
};

#endif // DELTAOPERATOR_H

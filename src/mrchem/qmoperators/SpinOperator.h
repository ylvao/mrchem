#ifndef SPINOPERATOR_H
#define SPINOPERATOR_H

#include "QMOperator.h"

class SpinOperator : public QMOperator {
public:
    SpinOperator() { }
    virtual ~SpinOperator() { }

    virtual void setup(double prec) { NOT_IMPLEMENTED_ABORT; }
    virtual void clear() { NOT_IMPLEMENTED_ABORT; }

    virtual Orbital* operator()(Orbital &orb_p) { NOT_IMPLEMENTED_ABORT; }
    virtual Orbital* adjoint(Orbital &orb_p) { NOT_IMPLEMENTED_ABORT; }

    using QMOperator::operator();
    using QMOperator::adjoint;
};

#endif // SPINOPERATOR_H

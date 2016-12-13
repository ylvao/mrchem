#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "QMOperator.h"
#include "ComplexFunction.h"
#include "OrbitalMultiplier.h"

class Potential : public ComplexFunction<3>, public QMOperator {
public:
    Potential();
    Potential(const Potential &pot) : QMOperator(pot) { NOT_IMPLEMENTED_ABORT; }
    Potential &operator=(const Potential &pot) { NOT_IMPLEMENTED_ABORT; }
    virtual ~Potential() { }

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

#endif // POTENTIAL_H

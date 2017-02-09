#ifndef QMPOTENTIAL_H
#define QMPOTENTIAL_H

#include "QMOperator.h"
#include "ComplexFunction.h"

class QMPotential : public ComplexFunction<3>, public QMOperator {
public:
    QMPotential();
    virtual ~QMPotential();

    virtual Orbital* operator() (Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    void calcRealPart(Orbital &Vphi, Orbital &phi);
    void calcImagPart(Orbital &Vphi, Orbital &phi, bool adjoint);
};

#endif // POTENTIAL_H

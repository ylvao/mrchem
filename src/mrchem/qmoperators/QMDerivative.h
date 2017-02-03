#ifndef QMDERIVATIVE_H
#define QMDERIVATIVE_H

#include "QMOperator.h"
#include "ABGVOperator.h"
#include "PHOperator.h"

class QMDerivative : public QMOperator {
public:
    QMDerivative(int dir, DerivativeOperator<3> &d);
    virtual ~QMDerivative() { }

    virtual void setup(double prec) { QMOperator::setup(prec); }
    virtual void clear() { QMOperator::clear(); }

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

protected:
    const int apply_dir;
    DerivativeOperator<3> *diff_oper;
};

#endif // QMDERIVATIVE_H

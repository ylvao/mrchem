#ifndef QMDERIVATIVE_H
#define QMDERIVATIVE_H

#include "QMOperator.h"
#include "ABGVOperator.h"
#include "PHOperator.h"

class QMDerivative : public QMOperator {
public:
    QMDerivative(int dir);
    virtual ~QMDerivative() { }

    virtual void setup(double prec) { QMOperator::setup(prec); }
    virtual void clear() { QMOperator::clear(); }

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

protected:
    const int apply_dir;
    ABGVOperator<3> diff_oper;
    //PHOperator<3> diff_oper;
};

#endif // QMDERIVATIVE_H

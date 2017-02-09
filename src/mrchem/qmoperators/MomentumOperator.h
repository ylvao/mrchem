#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMDerivative.h"

class MomentumOperator : public QMDerivative {
public:
    MomentumOperator(int dir, DerivativeOperator<3> &D)
        : QMDerivative(dir, D) { }
    virtual ~MomentumOperator() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

protected:
    virtual bool isReal() const { return false; }
    virtual bool isImag() const { return true; }
};

#endif // MOMENTUMOPERATOR_H


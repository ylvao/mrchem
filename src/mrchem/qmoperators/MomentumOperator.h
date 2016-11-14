#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMOperator.h"
#include "DerivativeOperator.h"
#include "OperatorApplier.h"

class MomentumOperator : public QMOperator {
public:
    MomentumOperator(int dir, double build_prec = -1.0);
    virtual ~MomentumOperator() { }

    virtual void setup(double prec);
    virtual void clear();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    const int apply_dir;
    DerivativeOperator<3> derivative;
    OperatorApplier<3> apply;
};

#endif // MOMENTUMOPERATOR_H


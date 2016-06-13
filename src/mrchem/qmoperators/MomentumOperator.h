#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMOperator.h"

template<int D> class DerivativeOperator;

class MomentumOperator : public QMOperator {
public:
    MomentumOperator(int dir, DerivativeOperator<3> &d_oper);
    virtual ~MomentumOperator();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    const int apply_dir;
    DerivativeOperator<3> *D;
};

#endif // MOMENTUMOPERATOR_H


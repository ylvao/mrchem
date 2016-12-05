#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMOperator.h"
#include "ABGVOperator.h"
#include "PHOperator.h"

class MomentumOperator : public QMOperator {
public:
    MomentumOperator(int dir, double prec = -1.0);
    virtual ~MomentumOperator() { }

    virtual void setup(double prec) { QMOperator::setup(prec); }
    virtual void clear() { QMOperator::clear(); }

    virtual int printTreeSizes() const;

    virtual Orbital* operator()(Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    const int apply_dir;
    ABGVOperator<3> diff_oper;
    //PHOperator<3> diff_oper;
};

#endif // MOMENTUMOPERATOR_H


#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMOperator.h"
#include "DerivativeOperator.h"
#include "GridGenerator.h"

template<int D> class DerivativeOperator;
template<int D> class MultiResolutionAnalysis;

class MomentumOperator : public QMOperator {
public:
    MomentumOperator(double build_prec,
                     const MultiResolutionAnalysis<3> &mra,
                     int dir);
    virtual ~MomentumOperator();

    virtual void setup(double prec);
    virtual void clear();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    GridGenerator<3> grid;
    DerivativeOperator<3> derivative;
};

#endif // MOMENTUMOPERATOR_H


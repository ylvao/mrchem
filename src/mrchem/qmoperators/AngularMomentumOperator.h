#ifndef ANGULARMOMENTUMOPERATOR_H
#define ANGULARMOMENTUMOPERATOR_H

#include "AngularOperator.h"

class AngularMomentumOperator : public AngularOperator {
public:
    AngularMomentumOperator(int d, const double *o)
        : AngularOperator(d, 1.0/2.0, o) {
    }
    virtual ~AngularMomentumOperator() { }
};

#endif // ANGULARMOMENTUMOPERATOR_H

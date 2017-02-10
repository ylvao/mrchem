#ifndef ANGULARMOMENTUMOPERATOR_H
#define ANGULARMOMENTUMOPERATOR_H

#include "QMTensorOperator.h"
#include "PositionOperator.h"
#include "MomentumOperator.h"
#include "DerivativeOperator.h"

class AngularMomentumOperator : public RankOneTensorOperator<3> {
public:
    AngularMomentumOperator(DerivativeOperator<3> &D, const double *o = 0)
            : r_x(0, o), r_y(1, o), r_z(2, o),
              p_x(0, D), p_y(1, D), p_z(2, D) {
        initializeTensorOperator();
    }
    virtual ~AngularMomentumOperator() { }

protected:
    PositionOperator r_x;
    PositionOperator r_y;
    PositionOperator r_z;
    MomentumOperator p_x;
    MomentumOperator p_y;
    MomentumOperator p_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (r_y*p_z - r_z*p_y);
        h[1] = (r_z*p_x - r_x*p_z);
        h[2] = (r_x*p_y - r_y*p_x);
    }
};

#endif // ANGULARMOMENTUMOPERATOR_H

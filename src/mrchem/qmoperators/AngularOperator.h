#ifndef ANGULAROPERATOR_H
#define ANGULAROPERATOR_H

#include "QMOperator.h"
#include "DipoleOperator.h";
#include "MomentumOperator.h";

class AngularOperator : public QMOperator {
public:
    AngularOperator(int dir, double c, const double *o = 0);
    virtual ~AngularOperator() { }

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator()(Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    int apply_dir;
    double coef;
    PositionOperator r_x;
    PositionOperator r_y;
    PositionOperator r_z;
    MomentumOperator p_x;
    MomentumOperator p_y;
    MomentumOperator p_z;
};

#endif // ANGULAROPERATOR_H



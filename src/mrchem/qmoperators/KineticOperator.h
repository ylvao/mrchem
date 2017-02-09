#ifndef KINETICOPERATOR_H
#define KINETICOPERATOR_H

#include "QMTensorOperator.h"
#include "MomentumOperator.h"

class KineticOperator : public RankZeroTensorOperator {
public:
    KineticOperator(DerivativeOperator<3> &D)
            : p_x(0, D), p_y(1, D), p_z(2, D) {
        initializeTensorOperator();
    }
    virtual ~KineticOperator() { }

    virtual void setup(double prec) {
        this->p_x.setup(prec);
        this->p_y.setup(prec);
        this->p_z.setup(prec);
    }
    virtual void clear() {
        this->p_x.clear();
        this->p_y.clear();
        this->p_z.clear();
    }

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::adjoint;

protected:
    MomentumOperator p_x;
    MomentumOperator p_y;
    MomentumOperator p_z;

    void initializeTensorOperator() {
        RankZeroTensorOperator &h = *this;
        h = -0.5*(p_x*p_x + p_y*p_y + p_z*p_z);
    }
};

#endif // KINETICOPERATOR_H

#ifndef KINETICOPERATOR_H
#define KINETICOPERATOR_H

#include "QMTensorOperator.h"
#include "MomentumOperator.h"

class KineticOperator : public RankOneTensorOperator<3> {
public:
    KineticOperator() : p() { initializeTensorOperator(); }
    virtual ~KineticOperator() { }

    virtual void setup(double prec) { this->p.setup(prec); }
    virtual void clear() { this->p.clear(); }

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    MomentumOperator p;

    void initializeTensorOperator() {
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankOneTensorOperator<3> &h = *this;
        h[0] = -0.5*p_x*p_x;
        h[1] = -0.5*p_y*p_z;
        h[2] = -0.5*p_y*p_z;
    }
};

#endif // KINETICOPERATOR_H

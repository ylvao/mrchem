#ifndef KINETICOPERATOR_H
#define KINETICOPERATOR_H

#include "QMOperator.h"

template<int D> class DerivativeOperator;
class MomentumOperator;
class Orbital;

class KineticOperator : public QMOperator {
public:
    KineticOperator(DerivativeOperator<3> &d);
    virtual ~KineticOperator();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    MomentumOperator *p_x;
    MomentumOperator *p_y;
    MomentumOperator *p_z;
};

#endif // KINETICOPERATOR_H

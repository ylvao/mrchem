#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "QMOperator.h"
#include "MWAdder.h"
#include "MWMultiplier.h"

template<int D> class FunctionTree;
template<int D> class MultiResolutionAnalysis;

class Potential : public QMOperator {
public:
    Potential(const MultiResolutionAnalysis<3> &mra);
    virtual ~Potential();

    int getNNodes() const;
    virtual int printTreeSizes() const;

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &orb);
    virtual Orbital* adjoint(Orbital &orb);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    FunctionTree<3> *real;
    FunctionTree<3> *imag;
};

#endif // POTENTIAL_H

#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "QMOperator.h"

template<int D> class MWAdder;
template<int D> class MWMultiplier;
template<int D> class FunctionTree;

class Potential : public QMOperator {
public:
    Potential(MWAdder<3> &a,
              MWMultiplier<3> &m,
              FunctionTree<3> *re = 0,
              FunctionTree<3> *im = 0)
        : add(&a),
          mult(&m),
          real(re),
          imag(im) {
    }
    virtual ~Potential() { }

    int getNNodes() const;
    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb);
    virtual Orbital* adjoint(Orbital &orb);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    MWAdder<3> *add;
    MWMultiplier<3> *mult;
    FunctionTree<3> *real;
    FunctionTree<3> *imag;
};

#endif // POTENTIAL_H

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
    Potential &operator=(const Potential &pot) { NOT_IMPLEMENTED_ABORT; }
    virtual ~Potential();

    int getNNodes() const;
    virtual int printTreeSizes() const;

    virtual void setup(double prec);
    virtual void clear();

    bool hasReal() const { if (this->real == 0) return false; return true; }
    bool hasImag() const { if (this->imag == 0) return false; return true; }

    FunctionTree<3> &re() { return *this->real; }
    FunctionTree<3> &im() { return *this->imag; }

    virtual Orbital* operator() (Orbital &orb);
    virtual Orbital* adjoint(Orbital &orb);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    friend class XCPotential;
    friend class CoulombPotential;

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    FunctionTree<3> *real;
    FunctionTree<3> *imag;
};

#endif // POTENTIAL_H

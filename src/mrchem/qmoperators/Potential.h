#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "QMOperator.h"
#include "OrbitalMultiplier.h"

template<int D> class FunctionTree;

class Potential : public QMOperator {
public:
    Potential();
    Potential &operator=(const Potential &pot) { NOT_IMPLEMENTED_ABORT; }
    virtual ~Potential();

    int getNNodes() const;
    virtual int printTreeSizes() const;

    virtual void setup(double prec);
    virtual void clear();

    bool hasReal() const { if (this->real == 0) return false; return true; }
    bool hasImag() const { if (this->imag == 0) return false; return true; }

    void allocReal();
    void allocImag();

    void setReal(FunctionTree<3> *re) { this->real = re; }
    void setImag(FunctionTree<3> *im) { this->imag = im; }

    FunctionTree<3> &re() { return *this->real; }
    FunctionTree<3> &im() { return *this->imag; }

    virtual Orbital* operator() (Orbital &orb);
    virtual Orbital* adjoint(Orbital &orb);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    OrbitalMultiplier mult;
    FunctionTree<3> *real;
    FunctionTree<3> *imag;
};

#endif // POTENTIAL_H

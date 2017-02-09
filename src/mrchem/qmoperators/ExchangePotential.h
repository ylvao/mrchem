#ifndef EXCHANGEPOTENTIAL_H
#define EXCHANGEPOTENTIAL_H

#include "ExchangeOperator.h"

class ExchangePotential : public ExchangeOperator {
public:
    ExchangePotential(PoissonOperator &P, OrbitalVector &phi, double x_fac = 1.0);
    virtual ~ExchangePotential() { }

    virtual void rotate(Eigen::MatrixXd &U);

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &phi_p);
    virtual Orbital* adjoint(Orbital &phi_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    OrbitalVector exchange;  ///< Set to keep precomputed exchange orbitals

    Orbital* calcExchange(Orbital &phi_p);
    Orbital* testPreComputed(Orbital &phi_p);
    void calcInternalExchange();
    int calcInternal(int i);
    int calcInternal(int i, int j);
};

#endif // EXCHANGEPOTENTIAL_H

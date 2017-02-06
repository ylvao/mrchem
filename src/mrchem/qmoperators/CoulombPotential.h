#ifndef COULOMBPOTENTIAL_H
#define COULOMBPOTENTIAL_H

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(PoissonOperator &P, OrbitalVector &phi)
        : CoulombOperator(P, phi) { }
    virtual ~CoulombPotential() { }

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);
};

#endif // COULOMBPOTENTIAL_H

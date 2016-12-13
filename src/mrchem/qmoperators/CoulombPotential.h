#ifndef COULOMBPOTENTIAL_H
#define COULOMBPOTENTIAL_H

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(PoissonOperator &P, OrbitalVector &phi)
        : CoulombOperator(P, phi) { }
    virtual ~CoulombPotential() { }

    void setup(double prec);
    void clear();
};

#endif // COULOMBPOTENTIAL_H

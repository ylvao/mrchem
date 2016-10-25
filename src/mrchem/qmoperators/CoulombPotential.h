#ifndef COULOMBPOTENTIAL_H
#define COULOMBPOTENTIAL_H

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(double build_prec, OrbitalVector &phi)
        : CoulombOperator(build_prec, phi) { }
    virtual ~CoulombPotential() { }

    void setup(double prec);
    void clear();
};

#endif // COULOMBPOTENTIAL_H

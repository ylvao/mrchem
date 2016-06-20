#ifndef COULOMBPOTENTIAL_H
#define COULOMBPOTENTIAL_H

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(double build_prec,
                    const MultiResolutionAnalysis<3> &mra,
                    OrbitalVector &phi);

    void setup(double prec);
    void clear();
};

#endif // COULOMBPOTENTIAL_H

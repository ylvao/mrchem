#ifndef HARTREE_H
#define HARTREE_H

#include "FockOperator.h"

class Hartree : public FockOperator {
public:
    Hartree(const MultiResolutionAnalysis<3> &mra,
            KineticOperator &t,
            NuclearPotential &v,
            CoulombOperator &j)
        : FockOperator(mra, &t, &v, &j) { }
    virtual ~Hartree() { }
};

#endif // HARTREE_H

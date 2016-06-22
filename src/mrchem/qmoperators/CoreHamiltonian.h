#ifndef COREHAMILTONIAN_H
#define COREHAMILTONIAN_H

#include "FockOperator.h"

class CoreHamiltonian : public FockOperator {
public:
    CoreHamiltonian(const MultiResolutionAnalysis<3> &mra,
                    KineticOperator &t,
                    NuclearPotential &v)
        : FockOperator(mra, &t, &v) { }
    virtual ~CoreHamiltonian() { }
};

#endif // COREHAMILTONIAN_H

#ifndef COREHAMILTONIAN_H
#define COREHAMILTONIAN_H

#include "FockOperator.h"

class CoreHamiltonian : public FockOperator {
public:
    CoreHamiltonian(KineticOperator &t, NuclearPotential &v)
        : FockOperator(&t, &v) { }
    virtual ~CoreHamiltonian() { }
};

#endif // COREHAMILTONIAN_H

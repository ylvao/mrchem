#ifndef COULOMBOPERATOR_H
#define COULOMBOPERATOR_H

#include "QMPotential.h"
#include "Density.h"

class CoulombOperator : public QMPotential {
public:
    CoulombOperator(PoissonOperator &P, OrbitalVector &phi)
        : poisson(&P),
          orbitals(&phi),
          density(false) {
    }
    virtual ~CoulombOperator() { }

protected:
    PoissonOperator *poisson;   // Pointer to external object
    OrbitalVector *orbitals;    // Pointer to external object

    Density density;            // Density that defines the potential
};

#endif // COULOMBOPERATOR_H

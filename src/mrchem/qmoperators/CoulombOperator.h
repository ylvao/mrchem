#ifndef COULOMBOPERATOR_H
#define COULOMBOPERATOR_H

#include "TwoElectronOperator.h"
#include "Density.h"
#include "QMPotential.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class CoulombOperator : public TwoElectronOperator {
public:
    CoulombOperator(PoissonOperator &P, OrbitalVector &phi)
        : TwoElectronOperator(MRA->getMaxScale(), phi),
          poisson(&P),
          density(false),
          potential() {
    }
    virtual ~CoulombOperator() { }

protected:
    PoissonOperator *poisson;   // Pointer to external object

    Density density;            // Density that defines the potential
    QMPotential potential;      // The actual operator
};

#endif // COULOMBOPERATOR_H

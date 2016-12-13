#ifndef COULOMBOPERATOR_H
#define COULOMBOPERATOR_H

#include "QMOperator.h"
#include "Density.h"
#include "Potential.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class CoulombOperator : public QMOperator {
public:
    CoulombOperator(PoissonOperator &P, OrbitalVector &phi)
        : QMOperator(MRA->getMaxScale()),
          poisson(&P),
          orbitals(&phi),
          density(false),
          potential() { }
    virtual ~CoulombOperator() { }

    virtual Orbital* operator() (Orbital &phi_p) {
        if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
        Potential &V = this->potential;
        return V(phi_p);
    }
    virtual Orbital* adjoint(Orbital &phi_p) {
        if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
        Potential &V = this->potential;
        return V.adjoint(phi_p);
    }

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    PoissonOperator *poisson;   // Pointer to external object
    OrbitalVector *orbitals;    // Pointer to external object
    Density density;            // Density that defines the potential
    Potential potential;        // The actual operator
};

#endif // COULOMBOPERATOR_H

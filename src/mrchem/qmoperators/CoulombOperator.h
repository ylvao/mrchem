#ifndef COULOMBOPERATOR_H
#define COULOMBOPERATOR_H

#include "QMOperator.h"
#include "Density.h"
#include "Potential.h"
#include "PoissonOperator.h"
#include "DensityProjector.h"

class OrbitalVector;

class CoulombOperator : public QMOperator {
public:
    CoulombOperator(double build_prec, OrbitalVector &phi);
    virtual ~CoulombOperator();

    virtual void setup(double prec);
    virtual void clear();

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    using QMOperator::operator();
    using QMOperator::adjoint;

protected:
    PoissonOperator poisson;
    DensityProjector project;

    Density density;
    Potential potential;
    OrbitalVector *orbitals;
};

#endif // COULOMBOPERATOR_H

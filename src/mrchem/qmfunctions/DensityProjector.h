#ifndef DENSITYPROJECTOR_H
#define DENSITYPROJECTOR_H

#include "GridGenerator.h"
#include "MWAdder.h"
#include "MWMultiplier.h"

class OrbitalVector;
class Orbital;
class Density;

class DensityProjector {
public:
    DensityProjector(double prec = -1.0);
    virtual ~DensityProjector() { }

    void setPrecision(double prec);

    void operator()(Density &rho, Orbital &phi);
    void operator()(Density &rho, OrbitalVector &phi);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;
};

#endif // DENSITYPROJECTOR_H

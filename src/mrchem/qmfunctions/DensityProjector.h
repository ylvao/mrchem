#ifndef DENSITYPROJECTOR_H
#define DENSITYPROJECTOR_H

#include "GridGenerator.h"
#include "GridCleaner.h"
#include "MWAdder.h"
#include "MWMultiplier.h"

class OrbitalVector;
class Orbital;
class Density;

class DensityProjector {
public:
    DensityProjector(const MultiResolutionAnalysis<3> &mra, double pr = -1.0)
        : grid(mra),
          clean(mra, pr),
          add(mra, pr),
          mult(mra, pr) { }
    virtual ~DensityProjector() { }

    void setPrecision(double prec);

    void operator()(Density &rho, Orbital &phi);
    void operator()(Density &rho, OrbitalVector &phi);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;
    GridCleaner<3> clean;
};

#endif // DENSITYPROJECTOR_H

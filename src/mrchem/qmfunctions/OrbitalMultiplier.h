#ifndef ORBITALMULTIPLIER_H
#define ORBITALMULTIPLIER_H

#include "MWMultiplier.h"
#include "MWAdder.h"
#include "Orbital.h"

class OrbitalMultiplier {
public:
    OrbitalMultiplier(const MultiResolutionAnalysis<3> &mra, double pr = -1.0);
    virtual ~OrbitalMultiplier() { }

    void setPrecision(double prec);

    void operator()(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);
    void adjoint(Orbital &phi_ab, double c, Orbital &phi_a, Orbital &phi_b);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;

    FunctionTree<3> *calcRealPart(double c, Orbital &phi_a, Orbital &phi_b);
    FunctionTree<3> *calcImagPart(double c, Orbital &phi_a, Orbital &phi_b, bool adjoint);
};


#endif // ORBITALMULTIPLIER_H

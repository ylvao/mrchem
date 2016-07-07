#ifndef ORBITALPROJECTOR_H
#define ORBITALPROJECTOR_H

#include <string>

#include "MWProjector.h"

class OrbitalVector;
class OrbitalExp;
class Nuclei;

class OrbitalProjector {
public:
    OrbitalProjector(const MultiResolutionAnalysis<3> &mra, double prec)
        : grid(mra),
          project(mra, prec) { }
    virtual ~OrbitalProjector() { }

    OrbitalVector* operator()(const Nuclei &nucs);

    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo);
    void operator()(OrbitalVector &phi,
                    const std::string &bf,
                    const std::string &mo_a,
                    const std::string &mo_b);
protected:
    GridGenerator<3> grid;
    MWProjector<3> project;

    OrbitalExp* readOrbitalExpansion(const std::string &bf,
                                     const std::string &mo);
};

#endif // ORBITALPROJECTOR_H

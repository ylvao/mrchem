#ifndef INITIALGUESSPROJECTOR_H
#define INITIALGUESSPROJECTOR_H

#include <string>

#include "MWProjector.h"

class OrbitalVector;
class OrbitalExp;

class InitialGuessProjector {
public:
    InitialGuessProjector(MultiResolutionAnalysis<3> &mra, double prec)
        : G(mra),
          Q(mra, prec) { }
    virtual ~InitialGuessProjector() { }

//    void readOrbitals(const OrbitalVector &orbs);
//    void readVirtuals(const std::string &bf,
//                      const std::string &mo,
//                      int n_occ);
    void operator()(OrbitalVector &orbs,
                    const std::string &bf,
                    const std::string &mo);
    void operator()(OrbitalVector &orbs,
                    const std::string &bf,
                    const std::string &mo_a,
                    const std::string &mo_b);
protected:
    GridGenerator<3> G;
    MWProjector<3> Q;

    OrbitalExp* readOrbitalExpansion(const std::string &bf,
                                     const std::string &mo);
};

#endif // INITIALGUESSPROJECTOR_H

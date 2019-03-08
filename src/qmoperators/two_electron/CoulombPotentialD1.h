#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD1 final : public CoulombPotential {
public:
    CoulombPotentialD1(mrcpp::PoissonOperator *P, OrbitalVector *Phi);

private:
    void setupLocalDensity(double prec) override;
    void setupGlobalDensity(double prec) override;
};

} // namespace mrchem

#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD2 final : public CoulombPotential {
public:
    CoulombPotentialD2(std::shared_ptr<mrcpp::PoissonOperator> P,
                       OrbitalVector *Phi,
                       OrbitalVector *X,
                       OrbitalVector *Y);

private:
    OrbitalVector *orbitals_x; ///< Perturbed orbitals
    OrbitalVector *orbitals_y; ///< Perturbed orbitals

    void setupGlobalDensity(double prec) override;
    void setupLocalDensity(double prec) override;
};

} // namespace mrchem

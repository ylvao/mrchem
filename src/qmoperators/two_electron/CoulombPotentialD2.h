#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD2 final : public CoulombPotential {
public:
    CoulombPotentialD2(std::shared_ptr<mrcpp::PoissonOperator> P,
                       std::shared_ptr<OrbitalVector> Phi,
                       std::shared_ptr<OrbitalVector> X,
                       std::shared_ptr<OrbitalVector> Y,
                       bool mpi_share = false);

private:
    std::shared_ptr<OrbitalVector> orbitals_x; ///< Perturbed orbitals
    std::shared_ptr<OrbitalVector> orbitals_y; ///< Perturbed orbitals

    void setupGlobalDensity(double prec) override;
    void setupLocalDensity(double prec) override;
};

} // namespace mrchem

#pragma once

#include "XCPotential.h"
#include "XCPotentialD1.h"
#include "XCPotentialD2.h"
#include "qmoperators/RankZeroTensorOperator.h"

/** @class XCOperator
 *
 * @brief DFT Exchange-Correlation operator containing a single XCPotential
 *
 * This class is a simple TensorOperator realization of @class XCPotential.
 *
 */

namespace mrchem {

class XCOperator final : public RankZeroTensorOperator {
public:
    explicit XCOperator(std::unique_ptr<mrdft::MRDFT> &F,
                        std::shared_ptr<OrbitalVector> Phi = nullptr,
                        bool mpi_shared = false) {
        potential = std::make_shared<XCPotentialD1>(F, Phi, mpi_shared);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &XC = (*this);
        XC = potential;
        XC.name() = "V_xc";
    }
    XCOperator(std::unique_ptr<mrdft::MRDFT> &F,
               std::shared_ptr<OrbitalVector> Phi,
               std::shared_ptr<OrbitalVector> X,
               std::shared_ptr<OrbitalVector> Y,
               bool mpi_shared = false) {
        potential = std::make_shared<XCPotentialD2>(F, Phi, X, Y, mpi_shared);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &XC = (*this);
        XC = potential;
        XC.name() = "V_xc";
    }
    ~XCOperator() override = default;

    auto getEnergy() { return potential->getEnergy(); }
    auto &getDensity(DensityType spin, int pert_idx = 0) { return potential->getDensity(spin, pert_idx); }

private:
    std::shared_ptr<XCPotential> potential{nullptr};
};

} // namespace mrchem

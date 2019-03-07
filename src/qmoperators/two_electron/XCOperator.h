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
    XCOperator(std::shared_ptr<mrdft::XCFunctional> F, std::shared_ptr<OrbitalVector> Phi = nullptr)
            : potential(new XCPotentialD1(F, Phi)) {
        RankZeroTensorOperator &XC = (*this);
        XC = *this->potential;
    }
    XCOperator(std::shared_ptr<mrdft::XCFunctional> F,
               std::shared_ptr<OrbitalVector> Phi,
               std::shared_ptr<OrbitalVector> X,
               std::shared_ptr<OrbitalVector> Y)
            : potential(new XCPotentialD2(F, Phi, X, Y)) {
        RankZeroTensorOperator &XC = (*this);
        XC = *this->potential;
    }
    ~XCOperator() override = default;

    void setupDensity(double prec = -1.0) { this->potential->setupDensity(prec); }
    void setupPotential(double prec = -1.0) { this->potential->setupPotential(prec); }

    double getEnergy() { return this->potential->getEnergy(); }
    mrcpp::FunctionTree<3> &getDensity(int spin) { return this->potential->getDensity(spin); }
    std::shared_ptr<mrdft::XCFunctional> &getFunctional() { return this->potential->getFunctional(); }

private:
    std::shared_ptr<XCPotential> potential;
};

} // namespace mrchem

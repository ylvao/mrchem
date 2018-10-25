#pragma once

#include "qmoperators/RankZeroTensorOperator.h"
#include "XCPotential.h"
#include "XCPotentialD1.h"

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
    XCOperator(mrdft::XCFunctional *F, OrbitalVector *Phi = nullptr)
            : potential(nullptr) {
        this->potential = new XCPotentialD1(F, Phi);
        RankZeroTensorOperator &XC = (*this);
        XC = *this->potential;
    }

    double getEnergy() { return this->potential->getEnergy(); }
    mrcpp::FunctionTree<3> &getDensity(int spin) { return this->potential->getDensity(spin); }

protected:
    XCPotential * potential;
};

} //namespace mrchem

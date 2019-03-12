#pragma once

#include "AngularMomentumOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class H_B_dip final : public RankOneTensorOperator<3> {
public:
    H_B_dip(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : l(D, o) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = 0.5 * l[0];
        h[1] = 0.5 * l[1];
        h[2] = 0.5 * l[2];
    }

private:
    AngularMomentumOperator l;
};

} // namespace mrchem

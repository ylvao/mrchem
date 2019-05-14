#pragma once

#include "MomentumOperator.h"
#include "PositionOperator.h"

namespace mrchem {

class AngularMomentumOperator final : public RankOneTensorOperator<3> {
public:
    AngularMomentumOperator(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &o)
            : r(o)
            , p(D) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = (r[1] * p[2] - r[2] * p[1]);
        h[1] = (r[2] * p[0] - r[0] * p[2]);
        h[2] = (r[0] * p[1] - r[1] * p[0]);
        h[0].name() = "l[x]";
        h[1].name() = "l[y]";
        h[2].name() = "l[z]";
    }

private:
    PositionOperator r;
    MomentumOperator p;
};

} // namespace mrchem

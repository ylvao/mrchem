#pragma once

#include "AngularMomentumOperator.h"
#include "DistanceOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class H_M_pso final : public RankOneTensorOperator<3> {
public:
    H_M_pso(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const mrcpp::Coord<3> &k)
            : r_m3(3.0, k)
            , l(D, k) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = alpha_2 * r_m3 * l[0];
        h[1] = alpha_2 * r_m3 * l[1];
        h[2] = alpha_2 * r_m3 * l[2];
        h[0].name() = "h_M_pso[x]";
        h[1].name() = "h_M_pso[y]";
        h[2].name() = "h_M_pso[z]";
    }

private:
    DistanceOperator r_m3;
    AngularMomentumOperator l;
};

} // namespace mrchem

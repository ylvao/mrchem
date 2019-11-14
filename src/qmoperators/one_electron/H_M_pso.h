#pragma once

#include "MomentumOperator.h"
#include "NuclearGradientOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class H_M_pso final : public RankOneTensorOperator<3> {
public:
    H_M_pso(std::shared_ptr<mrcpp::DerivativeOperator<3>> D, const Nucleus &nuc, double c)
            : p(D)
            , r_m3(nuc, c) {
        const double alpha_2 = PHYSCONST::alpha * PHYSCONST::alpha;

        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = alpha_2 * (r_m3[1] * p[2] - r_m3[2] * p[1]);
        h[1] = alpha_2 * (r_m3[2] * p[0] - r_m3[0] * p[2]);
        h[2] = alpha_2 * (r_m3[0] * p[1] - r_m3[1] * p[0]);
        h[0].name() = "h_M_pso[x]";
        h[1].name() = "h_M_pso[y]";
        h[2].name() = "h_M_pso[z]";
    }

private:
    MomentumOperator p;
    NuclearGradientOperator r_m3;
};

} // namespace mrchem

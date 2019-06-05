#pragma once

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/DistanceOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/PositionOperator.h"

namespace mrchem {

class X_rm3 final : public RankOneTensorOperator<3> {
public:
    X_rm3(const Nuclei &nucs, const mrcpp::Coord<3> &R_k)
            : r_m3(3.0, R_k, 1.0e-3)
            , r(nucs[0].getCoord()) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_m3 * r[0];
        h[1] = r_m3 * r[1];
        h[2] = r_m3 * r[2];
        h[0].name() = "x/r^3";
        h[1].name() = "y/r^3";
        h[2].name() = "z/r^3";
    }

private:
    DistanceOperator r_m3;
    PositionOperator r;
};

} // namespace mrchem

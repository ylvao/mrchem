#pragma once

#include "PositionOperator.h"
#include "qmoperators/RankTwoTensorOperator.h"

namespace mrchem {

class H_BB_dia final : public RankTwoTensorOperator<3, 3> {
public:
    H_BB_dia(const mrcpp::Coord<3> &o)
            : r(o) {
        // Invoke operator= to assign *this operator
        RankTwoTensorOperator<3, 3> &h = (*this);
        h[0][0] = (1.0 / 4.0) * (r[1] * r[1] + r[2] * r[2]);
        h[0][1] = (-1.0 / 4.0) * (r[0] * r[1]);
        h[0][2] = (-1.0 / 4.0) * (r[0] * r[2]);
        h[1][0] = (-1.0 / 4.0) * (r[1] * r[0]);
        h[1][1] = (1.0 / 4.0) * (r[0] * r[0] + r[2] * r[2]);
        h[1][2] = (-1.0 / 4.0) * (r[1] * r[2]);
        h[2][0] = (-1.0 / 4.0) * (r[2] * r[0]);
        h[2][1] = (-1.0 / 4.0) * (r[2] * r[1]);
        h[2][2] = (1.0 / 4.0) * (r[0] * r[0] + r[1] * r[1]);
        h[0][0].name() = "h_BB_dia[x,x]";
        h[0][1].name() = "h_BB_dia[x,y]";
        h[0][2].name() = "h_BB_dia[x,z]";
        h[1][0].name() = "h_BB_dia[y,x]";
        h[1][1].name() = "h_BB_dia[y,y]";
        h[1][2].name() = "h_BB_dia[y,z]";
        h[2][0].name() = "h_BB_dia[z,x]";
        h[2][1].name() = "h_BB_dia[z,y]";
        h[2][2].name() = "h_BB_dia[z,z]";
    }

private:
    PositionOperator r;
};

} // namespace mrchem

#pragma once

#include "PositionOperator.h"
#include "qmoperators/RankTwoTensorOperator.h"

/** Class H_E_quad
 *
 * @brief Electric quadrupole operator
 *
 * Q_ij = -3/2*r_i*r_j + d_ij/2*r^2
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges.
 *
 */
namespace mrchem {

class H_E_quad final : public RankTwoTensorOperator<3, 3> {
public:
    H_E_quad(const mrcpp::Coord<3> &o)
            : r(o) {
        // Invoke operator= to assign *this operator
        RankTwoTensorOperator &h = (*this);
        h[0][0] = -1.0 * r[0] * r[0] + 0.5 * r[1] * r[1] + 0.5 * r[2] * r[2];
        h[0][1] = -1.5 * r[0] * r[1];
        h[0][2] = -1.5 * r[0] * r[2];
        h[1][0] = -1.5 * r[1] * r[0];
        h[1][1] = -1.0 * r[1] * r[1] + 0.5 * r[0] * r[0] + 0.5 * r[2] * r[2];
        h[1][2] = -1.5 * r[1] * r[2];
        h[2][0] = -1.5 * r[2] * r[0];
        h[2][1] = -1.5 * r[2] * r[1];
        h[2][2] = -1.0 * r[2] * r[2] + 0.5 * r[0] * r[0] + 0.5 * r[1] * r[1];
        h[0][0].name() = "h_e_quad[x,x]";
        h[0][1].name() = "h_e_quad[x,y]";
        h[0][2].name() = "h_e_quad[x,z]";
        h[1][0].name() = "h_e_quad[y,x]";
        h[1][1].name() = "h_e_quad[y,y]";
        h[1][2].name() = "h_e_quad[y,z]";
        h[2][0].name() = "h_e_quad[z,x]";
        h[2][1].name() = "h_e_quad[z,y]";
        h[2][2].name() = "h_e_quad[z,z]";
    }

private:
    PositionOperator r;
};

} // namespace mrchem

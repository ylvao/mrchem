#pragma once

#include "PositionOperator.h"
#include "qmoperators/RankOneTensorOperator.h"

/** Class H_E_dip
 *
 * @brief Electric dipole operator
 *
 * mu_i = -r_i
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges.
 *
 */

namespace mrchem {

class H_E_dip final : public RankOneTensorOperator<3> {
public:
    H_E_dip(const mrcpp::Coord<3> &o)
            : r(o) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = -1.0 * r[0];
        h[1] = -1.0 * r[1];
        h[2] = -1.0 * r[2];
        h[0].name() = "h_E_dip[x]";
        h[1].name() = "h_E_dip[y]";
        h[2].name() = "h_E_dip[z]";
    }

protected:
    PositionOperator r;
};

} // namespace mrchem

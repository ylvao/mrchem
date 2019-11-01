#pragma once

#include "PositionOperator.h"
#include "chemistry/Nucleus.h"

/** Class H_E_dip
 *
 * @brief Electric dipole operator
 *
 * The operator is based on the PositionOperator. The sign is coherent
 * with nuclear and electronic charges. It impleents explicit trace
 * functions for the nuclear contributions.
 *
 */

namespace mrchem {

/** @brief constructor
 *
 * @param[in] o the origin wrt which the dipole moment is computed.
 *
 */
class H_E_dip final : public PositionOperator {
public:
    H_E_dip(const mrcpp::Coord<3> &o)
            : PositionOperator(o) {
        // Invoke operator= to assign *this operator
        RankOneTensorOperator<3> &h = (*this);
        h[0] = -1.0 * (*this)[0];
        h[1] = -1.0 * (*this)[1];
        h[2] = -1.0 * (*this)[2];
        h[0].name() = "h_E_dip[x]";
        h[1].name() = "h_E_dip[y]";
        h[2].name() = "h_E_dip[z]";
    }
};

} // namespace mrchem

/*
 *
 *
 *  \date Oct 15, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef SCALINGBASIS_H
#define SCALINGBASIS_H

#include "TelePrompter.h"
#include "constants.h"

class ScalingBasis {
public:
    ScalingBasis(int k, int t) : type(t), order(k) {
        if (this->order < 1) MSG_FATAL("Invalid scaling order");
    }
    virtual ~ScalingBasis() { }

    int getType() const { return this->type; }
    int getScalingOrder() const { return this->order; }
    int getQuadratureOrder() const { return this->order + 1; }

    bool operator==(const ScalingBasis &basis) const {
        if (this->type != basis.type) return false;
        if (this->order != basis.order) return false;
        return true;
    }
    bool operator!=(const ScalingBasis &basis) const {
        if (this->type != basis.type) return true;
        if (this->order != basis.order) return true;
        return false;
    }

    friend std::ostream& operator<<(std::ostream &o, const ScalingBasis &bas) {
        o << "*ScalingBasis:" << std::endl;
        o << "  order           = " << bas.getScalingOrder() << std::endl;
        if (bas.getType() == Legendre) {
            o << "  type            = Legendre";
        } else if (bas.getType() == Interpol) {
            o << "  type            = Interpolating";
        } else {
            o << "  type            = Unknown";
        }
        return o;
    }
protected:
    const int type;
    const int order;
};

#endif /* SCALINGBASIS_H */

/*
 *
 *
 *  \date Jul 26, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "QuadratureCache.h"

QuadratureCache::QuadratureCache() {
    A = 0.0;
    B = 1.0;
    intervals = 1;
}

QuadratureCache::~QuadratureCache() {
}

void QuadratureCache::load(int order) {
    SET_CACHE_LOCK();
    if (not hasId(order)) {
        GaussQuadrature *gp = new GaussQuadrature(order, A, B, intervals);
        int memo = 2 * order * sizeof(double);
        ObjectCache<GaussQuadrature>::load(order, gp, memo);
    }
    UNSET_CACHE_LOCK();
}

GaussQuadrature &QuadratureCache::get(int order) {
    if (not hasId(order)) {
        load(order);
    }
    return ObjectCache<GaussQuadrature>::get(order);
}

void QuadratureCache::setBounds(double a, double b) {
    if (fabs(A - a) < MachineZero and fabs(B - b) < MachineZero) {
        return;
    }
    if (a >= b) {
        MSG_ERROR("Invalid Gauss interval, a > b.");
    }
    A = a;
    B = b;
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) {
            ObjectCache<GaussQuadrature>::get(i).setBounds(a, b);
        }
    }
}

void QuadratureCache::setIntervals(int ivals) {
    if (ivals == intervals) {
        return;
    }
    if (intervals < 1) {
        MSG_ERROR("Invalid number of intervals, intervals < 1");
    }
    for (int i = 0; i < getNObjs(); i++) {
        if (hasId(i)) {
            ObjectCache<GaussQuadrature>::get(i).setIntervals(ivals);
        }
    }

}

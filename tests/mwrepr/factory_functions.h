#ifndef FACTORY_FUNCTIONS_H
#define FACTORY_FUNCTIONS_H

#include "BoundingBox.h"
#include "NodeIndex.h"
#include "MRTree.h"

template<class T> void finalize(T **obj) {
    if (obj == 0) MSG_FATAL("Invalid argument");
    if (*obj == 0) MSG_FATAL("Invalid argument");
    delete *obj;
    *obj = 0;
}

template<int D> void initialize(NodeIndex<D> **idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");
    if (*idx != 0) MSG_FATAL("Invalid argument");
    int scale = 1;
    int rank = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = d-1;
    }
    *idx = new NodeIndex<D>(scale, l, rank);
}

template<int D> void testInitial(const NodeIndex<D> *idx) {
    if (idx == 0) MSG_FATAL("Invalid argument");

    const int scale = 1;
    REQUIRE( scale == idx->getScale() );

    const int rank = 1;
    REQUIRE( rank == idx->getRankId() );

    for (int d = 0; d < D; d++) {
        const int l = d-1;
        REQUIRE( l == idx->getTranslation(d) );
        REQUIRE( l == idx->getTranslation()[d] );
    }
}

template<int D> void initialize(BoundingBox<D> **box) {
    if (box == 0) MSG_FATAL("Invalid argument");
    if (*box != 0) MSG_FATAL("Invalid argument");

    int nb[D];
    for (int d = 0; d < D; d++) {
        nb[d] = D-d;
    }
    NodeIndex<D> *nIdx = 0;
    initialize<D>(&nIdx);

    *box = new BoundingBox<D>(*nIdx, nb);
}

template<int D> void testInitial(const BoundingBox<D> *box) {
    if (box == 0) MSG_FATAL("Invalid argument");

    const NodeIndex<D> &cIdx = box->getCornerIndex();
    testInitial<D>(&cIdx);

    REQUIRE( box->getUnitLength() > 0.0 );

    for (int d = 0; d < D; d++) {
        REQUIRE( box->getBoxLength(d) > 0.0 );
        REQUIRE( box->getBoxLengths()[d] > 0.0 );

        REQUIRE( box->getLowerBound(d) < box->getUpperBound(d) );
        REQUIRE( box->getLowerBounds()[d] < box->getUpperBounds()[d] );
    }

    int tot_boxes = 1;
    for (int d = 0; d < D; d++) {
        const int nb = D-d;
        REQUIRE( nb == box->getNBoxes(d) );
        tot_boxes *= box->getNBoxes(d);
    }
    REQUIRE( tot_boxes == box->getNBoxes() );
}

#endif //FACTORY_FUNCTIONS_H

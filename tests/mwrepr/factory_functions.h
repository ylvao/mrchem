#ifndef FACTORY_FUNCTIONS_H
#define FACTORY_FUNCTIONS_H

#include "NodeBox.h"
#include "NodeIndex.h"
#include "MultiResolutionAnalysis.h"
#include "InterpolatingBasis.h"
#include "MWTree.h"

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
    initialize(&nIdx);

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
        REQUIRE( nb == box->size(d) );
        tot_boxes *= box->size(d);
    }
    REQUIRE( tot_boxes == box->size() );
}

template<int D> void initialize(NodeBox<D> **box) {
    if (box == 0) MSG_FATAL("Invalid argument");
    if (*box != 0) MSG_FATAL("Invalid argument");

    int nb[D];
    for (int d = 0; d < D; d++) {
        nb[d] = D-d;
    }
    NodeIndex<D> *nIdx = 0;
    initialize(&nIdx);

    *box = new NodeBox<D>(*nIdx, nb);
}

template<int D> void testInitial(const NodeBox<D> *box) {
    if (box == 0) MSG_FATAL("Invalid argument");

    REQUIRE( box->getNOccupied() == 0 );
}

template<int D> void initialize(MRTree<D> **tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");
    if (*tree != 0) MSG_FATAL("Invalid argument");

    BoundingBox<D> *world = 0;
    initialize(&world);

    *tree = new MRTree<D>(*world);
}

template<int D> void testInitial(MRTree<D> *tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");

    REQUIRE( tree->getDepth() == 1 );
    REQUIRE( tree->getNNodes() == 0 );
    REQUIRE( tree->getNEndNodes() == 0 );
    REQUIRE( tree->getNGenNodes() == 0 );
    REQUIRE( tree->getNAllocGenNodes() == 0 );
}

template<int D> void initialize(MWTree<D> **tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");
    if (*tree != 0) MSG_FATAL("Invalid argument");

    BoundingBox<D> *world = 0;
    initialize(&world);

    int k = 5;
    InterpolatingBasis basis(k);

    MultiResolutionAnalysis<D> mra(*world, basis);

    *tree = new MWTree<D>(mra);
}

template<int D> void testInitial(const MWTree<D> *tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");

    REQUIRE( tree->getSquareNorm() == Approx(0.0) );
    REQUIRE( tree->getOrder() == 5 );
}

#endif //FACTORY_FUNCTIONS_H

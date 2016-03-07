#ifndef FACTORY_FUNCTIONS_H
#define FACTORY_FUNCTIONS_H

#include "BoundingBox.h"
#include "NodeIndex.h"
#include "MultiResolutionAnalysis.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"
#include "FunctionTree.h"
#include "GridGenerator.h"

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
        nb[d] = d + 1;
    }
    NodeIndex<D> *nIdx = 0;
    initialize(&nIdx);

    *box = new BoundingBox<D>(*nIdx, nb);
    finalize(&nIdx);
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
        const int nb = d + 1;
        REQUIRE( box->size(d) == nb );
        tot_boxes *= box->size(d);
    }
    REQUIRE( box->size() == tot_boxes );
}

template<int D> void initialize(FunctionTree<D> **tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");
    if (*tree != 0) MSG_FATAL("Invalid argument");

    BoundingBox<D> *world = 0;
    initialize(&world);

    int k = 5;
    InterpolatingBasis basis(k);

    MultiResolutionAnalysis<D> mra(*world, basis);
    GridGenerator<D> G(mra);
    *tree = G();

    finalize(&world);
}

template<int D> void testInitial(FunctionTree<D> *tree) {
    if (tree == 0) MSG_FATAL("Invalid argument");

    int k = 5;
    int tot_nodes = 1;
    for (int d = 0; d < D; d++) {
        tot_nodes *= D-d;
    }

    REQUIRE( tree->getSquareNorm() == Approx(-1.0) );
    REQUIRE( tree->getOrder() == k );
    REQUIRE( tree->getDepth() == 1 );
    REQUIRE( tree->getNNodes() == tot_nodes );
    REQUIRE( tree->getNEndNodes() == tot_nodes );
    REQUIRE( tree->getNGenNodes() == 0 );
    REQUIRE( tree->getNAllocGenNodes() == 0 );
}

#endif //FACTORY_FUNCTIONS_H

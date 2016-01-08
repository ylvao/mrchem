#include "catch.hpp"

#include "factory_functions.h"

namespace node_box {

template<int D> void testConstructors();

TEST_CASE("NodeBox constructor", "[node_box_constructor], [node_box], [boxes]") {
    SECTION("1D") {
        testConstructors<1>();
    }
    SECTION("2D") {
        testConstructors<2>();
    }
    SECTION("3D") {
        testConstructors<3>();
    }
}

template<int D> void testConstructors() {
    int nb[D];
    int tot_boxes = 1;
    for (int d = 0; d < D; d++) {
        nb[d] = D+d;
        tot_boxes *= nb[d];
    }
    NodeIndex<D> *nIdx = 0;
    initialize(&nIdx);

    NodeBox<D> box(*nIdx, nb);
    finalize(&nIdx);

    SECTION("Constructor") {
        REQUIRE( box.size() == tot_boxes);
        REQUIRE( box.getNOccupied() == 0 );
    }

    SECTION("Copy constructor") {
        NodeBox<D> box_copy(box);
        REQUIRE( box_copy.size() == tot_boxes);
        REQUIRE( box_copy.getNOccupied() == 0 );
    }

    SECTION("Base class copy constructor") {
        const BoundingBox<D> &b_box = static_cast<const BoundingBox<D> &>(box);
        NodeBox<D> box_copy(b_box);
        REQUIRE( box_copy == b_box );
    }
}

} // namespace

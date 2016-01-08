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
    NodeBox<D> *box = 0;
    initialize<D>(&box);

    SECTION("Constructor") {
        testInitial<D>(box);
    }

    SECTION("Copy constructor") {
        NodeBox<D> *box_copy = new NodeBox<D>(*box);
        testInitial<D>(box_copy);
        finalize(&box_copy);
    }

    SECTION("Base class copy constructor") {
        const BoundingBox<D> *b_box = static_cast<const BoundingBox<D> *>(box);
        NodeBox<D> *box_copy = new NodeBox<D>(*b_box);
        testInitial<D>(box_copy);
        finalize(&box_copy);
    }

    finalize(&box);
}

} // namespace

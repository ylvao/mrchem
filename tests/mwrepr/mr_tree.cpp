#include "catch.hpp"

#include "factory_functions.h"

namespace mr_tree {

template<int D> void testConstructors();

TEST_CASE("MRTree constructors", "[mr_tree_constructor], [mr_tree], [trees]") {
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
    BoundingBox<D> *world = 0;
    initialize(&world);

    MRTree<D> tree(*world);
    finalize(&world);

    SECTION("Constructor") {
        REQUIRE( tree.getDepth() == 1 );
        REQUIRE( tree.getNNodes() == 0 );
        REQUIRE( tree.getNEndNodes() == 0 );
        REQUIRE( tree.getNGenNodes() == 0 );
        REQUIRE( tree.getNAllocGenNodes() == 0 );
    }

    SECTION("Copy constructor") {
        MRTree<D> tree_copy(tree);
        REQUIRE( tree_copy.getDepth() == 1 );
        REQUIRE( tree_copy.getNNodes() == 0 );
        REQUIRE( tree_copy.getNEndNodes() == 0 );
        REQUIRE( tree_copy.getNGenNodes() == 0 );
        REQUIRE( tree_copy.getNAllocGenNodes() == 0 );
    }
}

} // namespace

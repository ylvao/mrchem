#include "catch.hpp"

#include "factory_functions.h"

template<int D> void testConstructors();

TEST_CASE("MRTree constructors", "[mr_tree_constructor], [mr_tree], [trees]") {
    SECTION("1D") {
        REQUIRE( false );
    }
    SECTION("2D") {
        REQUIRE( false );
    }
    SECTION("3D") {
        REQUIRE( false );
    }
}

//template<int D> void testConstructors() {
//    MRTree<D> *tree = 0;
//    initialize(&tree);

//    SECTION("Constructor") {
//        testInitial(tree);
//    }

//    SECTION("Copy constructor") {
//        MRTree<D> *tree_copy = new MRTree<D>(*tree);
//        testInitial(tree_copy);
//        finalize(&tree_copy);
//    }

//    SECTION("Default constructor") {
//        MRTree<D> *tree_copy = new MRTree<D>();
//        REQUIRE( false );
//        SECTION("Assignment operator") {
//            *tree_copy = *tree;
//            REQUIRE( false );
//        }
//        finalize(&tree);
//    }

//    finalize(&tree);
//}

//template<int D> void initialize(MRTree<D> **tree) {
//    if (tree == 0) MSG_FATAL("Invalid argument");
//    if (*tree != 0) MSG_FATAL("Invalid argument");

//    BoundingBox<D> *world = 0;
//    initialize<D>(&world);

//    *tree = new MRTree<D>(*world);
//}

//template<int D> void testInitial(MRTree<D> *tree) {
//    if (tree == 0) MSG_FATAL("Invalid argument");

//    REQUIRE( tree->getDepth() == 0 );
//    REQUIRE( tree->getNNodes() == 0 );
//    REQUIRE( tree->getNRootNodes() == 0 );
//    REQUIRE( tree->getNEndNodes() == 0 );
//    REQUIRE( tree->getNGenNodes() == 0 );
//    REQUIRE( tree->getNAllocGenNodes() == 0 );
//}


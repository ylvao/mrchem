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
    MRTree<D> *tree = 0;
    initialize(&tree);

    SECTION("Constructor") {
        testInitial(tree);
    }

    SECTION("Copy constructor") {
        MRTree<D> *tree_copy = new MRTree<D>(*tree);
        testInitial(tree_copy);
        finalize(&tree_copy);
    }

    finalize(&tree);
}

} // namespace

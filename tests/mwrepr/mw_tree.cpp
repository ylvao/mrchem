#include "catch.hpp"

#include "factory_functions.h"

namespace mw_tree {

template<int D> void testConstructors();

TEST_CASE("MWTree constructors", "[mw_tree_constructor], [mw_tree], [trees]") {
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
    MWTree<D> *tree = 0;
    initialize(&tree);

    SECTION("Constructor") {
        testInitial(tree);
    }

    SECTION("Copy constructor") {
        MWTree<D> *tree_copy = new MWTree<D>(*tree);
        testInitial(tree_copy);
        finalize(&tree_copy);
    }

    finalize(&tree);
}

} // namespace

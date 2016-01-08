#include "catch.hpp"

#include "factory_functions.h"

namespace function_tree {

template<int D> void testConstructors();

TEST_CASE("FunctionTree constructor", "[function_tree_constructor], [function_tree], [trees]") {
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
    FunctionTree<D> *tree = 0;
    initialize(&tree);

    SECTION("Constructor") {
        testInitial(tree);
    }

    SECTION("Copy constructor") {
        FunctionTree<D> *tree_copy = new FunctionTree<D>(*tree);
        testInitial(tree_copy);
        finalize(&tree_copy);
    }

    SECTION("Base class copy constructor") {
        const MWTree<D> *mw_tree = static_cast<const MWTree<D> *>(tree);
        FunctionTree<D> *tree_copy = new FunctionTree<D>(*mw_tree);
        testInitial(tree_copy);
        finalize(&tree_copy);
    }

    finalize(&tree);
}

} // namespace

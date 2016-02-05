#include "catch.hpp"

#include "factory_functions.h"

namespace function_tree {

template<int D> void testConstructors();

TEST_CASE("FunctionTree: Constructor", "[function_tree_constructor], [function_tree], [trees]") {
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

SCENARIO("FunctionTree: Zero function", "[function_tree_zero], [function_tree], [trees]") {
    double r[3] = {-0.2, 0.6, 0.76};
    GIVEN("a zero function in 1D") {
        FunctionTree<1> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("its value in an arbitrary point is zero") {
            REQUIRE( tree->evalf(r) == Approx(0.0) );
        }
        THEN("it integrates to zero") {
            REQUIRE( tree->integrate() == Approx(0.0) );
        }
        THEN("the dot product with itself is zero") {
            REQUIRE( tree->dot(*tree) == Approx(0.0) );
        }
        finalize(&tree);
    }
    GIVEN("a zero function in 2D") {
        FunctionTree<2> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("its value in an arbitrary point is zero") {
            REQUIRE( tree->evalf(r) == Approx(0.0) );
        }
        THEN("it integrates to zero") {
            REQUIRE( tree->integrate() == Approx(0.0) );
        }
        THEN("the dot product with itself is zero") {
            REQUIRE( tree->dot(*tree) == Approx(0.0) );
        }
        finalize(&tree);
    }
    GIVEN("a zero function in 3D") {
        FunctionTree<3> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("its value in an arbitrary point is zero") {
            REQUIRE( tree->evalf(r) == Approx(0.0) );
        }
        THEN("it integrates to zero") {
            REQUIRE( tree->integrate() == Approx(0.0) );
        }
        THEN("the dot product with itself is zero") {
            REQUIRE( tree->dot(*tree) == Approx(0.0) );
        }
        finalize(&tree);
    }
}

} // namespace

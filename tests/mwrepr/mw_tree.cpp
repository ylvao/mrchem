#include "catch.hpp"

#include "factory_functions.h"
#include "MultiResolutionAnalysis.h"
#include "InterpolatingBasis.h"
#include "LegendreBasis.h"

namespace mw_tree {

template<int D> void testConstructors(const ScalingBasis &basis);

TEST_CASE("MWTree constructors", "[mw_tree_constructor], [mw_tree], [trees]") {
    const int k = 5;
    SECTION("Interpolating 1D") {
        InterpolatingBasis basis(k);
        testConstructors<1>(basis);
    }
    SECTION("Interpolating 2D") {
        InterpolatingBasis basis(k);
        testConstructors<2>(basis);
    }
    SECTION("Interpolating 3D") {
        InterpolatingBasis basis(k);
        testConstructors<3>(basis);
    }
    SECTION("Legendre 1D") {
        LegendreBasis basis(k);
        testConstructors<1>(basis);
    }
    SECTION("Legendre 2D") {
        LegendreBasis basis(k);
        testConstructors<2>(basis);
    }
    SECTION("Legendre 3D") {
        LegendreBasis basis(k);
        testConstructors<3>(basis);
    }
}

template<int D> void testConstructors(const ScalingBasis &basis) {
    BoundingBox<D> *world = 0;
    initialize(&world);

    MultiResolutionAnalysis<D> mra(*world, basis);
    finalize(&world);

    MWTree<D> tree(mra);

    SECTION("Constructor") {
        REQUIRE( tree.getSquareNorm() == Approx(-1.0) );
        REQUIRE( tree.getOrder() == 5 );
        REQUIRE( tree.getDepth() == 1 );
        REQUIRE( tree.getNNodes() == 0 );
        REQUIRE( tree.getNEndNodes() == 0 );
        REQUIRE( tree.getNGenNodes() == 0 );
        REQUIRE( tree.getNAllocGenNodes() == 0 );
    }

    SECTION("Copy constructor") {
        MWTree<D> tree_copy(tree);
        REQUIRE( tree_copy.getSquareNorm() == Approx(-1.0) );
        REQUIRE( tree_copy.getOrder() == 5 );
        REQUIRE( tree_copy.getDepth() == 1 );
        REQUIRE( tree_copy.getNNodes() == 0 );
        REQUIRE( tree_copy.getNEndNodes() == 0 );
        REQUIRE( tree_copy.getNGenNodes() == 0 );
        REQUIRE( tree_copy.getNAllocGenNodes() == 0 );
    }
}

} // namespace

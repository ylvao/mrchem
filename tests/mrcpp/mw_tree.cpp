#include "catch.hpp"

#include "factory_functions.h"

namespace mw_tree {

//template<int D> void testConstructors(const ScalingBasis &basis);
template<int D> void testNodeFetchers();

//TEST_CASE("MWTree: Constructors", "[mw_tree_constructor], [mw_tree], [trees]") {
//    const int k = 5;
//    SECTION("Interpolating 1D") {
//        InterpolatingBasis basis(k);
//        testConstructors<1>(basis);
//    }
//    SECTION("Interpolating 2D") {
//        InterpolatingBasis basis(k);
//        testConstructors<2>(basis);
//    }
//    SECTION("Interpolating 3D") {
//        InterpolatingBasis basis(k);
//        testConstructors<3>(basis);
//    }
//    SECTION("Legendre 1D") {
//        LegendreBasis basis(k);
//        testConstructors<1>(basis);
//    }
//    SECTION("Legendre 2D") {
//        LegendreBasis basis(k);
//        testConstructors<2>(basis);
//    }
//    SECTION("Legendre 3D") {
//        LegendreBasis basis(k);
//        testConstructors<3>(basis);
//    }
//}

//template<int D> void testConstructors(const ScalingBasis &basis) {
//    BoundingBox<D> *world = 0;
//    initialize(&world);

//    MultiResolutionAnalysis<D> mra(*world, basis);
//    finalize(&world);

//    MWTree<D> tree(mra);

//    SECTION("Constructor") {
//        REQUIRE( tree.getSquareNorm() == Approx(-1.0) );
//        REQUIRE( tree.getOrder() == 5 );
//        REQUIRE( tree.getDepth() == 1 );
//        REQUIRE( tree.getNNodes() == 0 );
//        REQUIRE( tree.getNEndNodes() == 0 );
//        REQUIRE( tree.getNGenNodes() == 0 );
//        REQUIRE( tree.getNAllocGenNodes() == 0 );
//    }

//    SECTION("Copy constructor") {
//        MWTree<D> tree_copy(tree);
//        REQUIRE( tree_copy.getSquareNorm() == Approx(-1.0) );
//        REQUIRE( tree_copy.getOrder() == 5 );
//        REQUIRE( tree_copy.getDepth() == 1 );
//        REQUIRE( tree_copy.getNNodes() == 0 );
//        REQUIRE( tree_copy.getNEndNodes() == 0 );
//        REQUIRE( tree_copy.getNGenNodes() == 0 );
//        REQUIRE( tree_copy.getNAllocGenNodes() == 0 );
//    }
//}


TEST_CASE("MWTree: Fetching nodes", "[mw_tree_fetch], [mw_tree], [trees]") {
    SECTION("1D") {
        testNodeFetchers<1>();
    }
    SECTION("2D") {
        testNodeFetchers<2>();
    }
    SECTION("3D") {
        testNodeFetchers<3>();
    }
}

template<int D> void testNodeFetchers() {
    const double r[3] = {-0.3, 0.6, 1.9};

    int cIdx = 1 << (D - 1);
    NodeIndex<D> *root = 0;
    initialize(&root);
    const NodeIndex<D> idx_0(*root);
    const NodeIndex<D> idx_1(idx_0, cIdx);
    const NodeIndex<D> idx_2(idx_1, cIdx);
    finalize(&root);

    FunctionTree<D> *tree = 0;
    initialize(&tree);
    tree->setZero();

    const FunctionTree<D> *const_tree = const_cast<const FunctionTree<D> *>(tree);

    // Fetch by NodeIndex
    SECTION("Find node by NodeIndex: existing node") {
        MWNode<D> *node = tree->findNode(idx_0);
        REQUIRE( node != 0 );
        REQUIRE( node->getNodeIndex() == idx_0 );
        REQUIRE( node->isAllocated() );
        REQUIRE( node->hasCoefs() );
        REQUIRE_FALSE( node->isGenNode() );
    }
    SECTION("Find const node by NodeIndex: existing node") {
        const MWNode<D> *node = const_tree->findNode(idx_0);
        REQUIRE( node != 0 );
        REQUIRE( node->getNodeIndex() == idx_0 );
        REQUIRE( node->isAllocated() );
        REQUIRE( node->hasCoefs() );
        REQUIRE_FALSE( node->isGenNode() );
    }
    SECTION("Find node by NodeIndex: non-existing node") {
        MWNode<D> *node = tree->findNode(idx_1);
        REQUIRE( node == 0 );
    }
    SECTION("Find const node by NodeIndex: non-existing node") {
        const MWNode<D> *node = const_tree->findNode(idx_2);
        REQUIRE( node == 0 );
    }
    SECTION("Get node by NodeIndex: existing node") {
        MWNode<D> &node = tree->getNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE_FALSE( node.isGenNode() );
    }
    SECTION("Get node by NodeIndex: non-existing node") {
        MWNode<D> &node = tree->getNode(idx_2);
        REQUIRE( node.getNodeIndex() == idx_2 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.isGenNode() );
    }
    SECTION("Get node or end node by NodeIndex: existing node") {
        MWNode<D> &node = tree->getNodeOrEndNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE_FALSE( node.isGenNode() );
    }
    SECTION("Get const node or end node by NodeIndex: existing node") {
        const MWNode<D> &node = const_tree->getNodeOrEndNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE_FALSE( node.isGenNode() );
    }
    SECTION("Get node or end node by NodeIndex: non-existing node") {
        MWNode<D> &node = tree->getNodeOrEndNode(idx_2);
        REQUIRE( node.getNodeIndex() != idx_2 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.isEndNode() );
        REQUIRE_FALSE( node.isGenNode() );
    }
    SECTION("Get const node or end node by NodeIndex: non-existing node") {
        const MWNode<D> &node = const_tree->getNodeOrEndNode(idx_1);
        REQUIRE( node.getNodeIndex() != idx_2 );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.isEndNode() );
        REQUIRE_FALSE( node.isGenNode() );
    }

    // Fetch by coordinate
    SECTION("Get node by coord: existing node") {
        int depth = 0;
        MWNode<D> &node = tree->getNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node by coord: non-existing node") {
        int depth = 3;
        MWNode<D> &node = tree->getNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node or end node by coord: existing node") {
        int depth = 0;
        MWNode<D> &node = tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get const node or end node by coord: existing node") {
        int depth = 0;
        const MWNode<D> &node = const_tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node or end node by coord: non-existing node") {
        int depth = 3;
        MWNode<D> &node = tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isEndNode() );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() != depth );
    }
    SECTION("Get const node or end node by coord: non-existing node") {
        int depth = 3;
        const MWNode<D> &node = const_tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isEndNode() );
        REQUIRE( node.isAllocated() );
        REQUIRE( node.hasCoefs() );
        REQUIRE( node.getDepth() != depth );
    }

    finalize(&tree);
}

SCENARIO("MWTree: Generating nodes", "[mw_tree_generating], [mw_tree], [trees]") {
    const double r[3] = {-0.3, 0.6, 1.9};
    const int depth = 3;
    GIVEN("a default function in 1D") {
        FunctionTree<1> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("there are no GenNodes") {
            REQUIRE( tree->getNGenNodes() == 0 );
            REQUIRE( tree->getNAllocGenNodes() == 0 );
        }
        WHEN("a non-existing node is fetched") {
            MWNode<1> &node = tree->getNode(r, depth);
            THEN("there will be allocated GenNodes") {
                REQUIRE( tree->getNGenNodes() > 0 );
                REQUIRE( tree->getNAllocGenNodes() > 0 );
            }
            AND_WHEN("the GenNodes are cleared") {
                tree->clearGenerated();
                THEN("there will be un-allocated GenNodes") {
                    REQUIRE( tree->getNGenNodes() > 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
            AND_WHEN("the GenNodes are deleted") {
                tree->deleteGenerated();
                THEN("there will be no GenNodes") {
                    REQUIRE( tree->getNGenNodes() == 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
        }
        finalize(&tree);
    }
    GIVEN("a default function in 2D") {
        FunctionTree<2> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("there are no GenNodes") {
            REQUIRE( tree->getNGenNodes() == 0 );
            REQUIRE( tree->getNAllocGenNodes() == 0 );
        }
        WHEN("a non-existing node is fetched") {
            MWNode<2> &node = tree->getNode(r, depth);
            THEN("there will be allocated GenNodes") {
                REQUIRE( tree->getNGenNodes() > 0 );
                REQUIRE( tree->getNAllocGenNodes() > 0 );
            }
            AND_WHEN("the GenNodes are cleared") {
                tree->clearGenerated();
                THEN("there will be un-allocated GenNodes") {
                    REQUIRE( tree->getNGenNodes() > 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
            AND_WHEN("the GenNodes are deleted") {
                tree->deleteGenerated();
                THEN("there will be no GenNodes") {
                    REQUIRE( tree->getNGenNodes() == 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
        }
        finalize(&tree);
    }
    GIVEN("a default function in 3D") {
        FunctionTree<3> *tree = 0;
        initialize(&tree);
        tree->setZero();
        THEN("there are no GenNodes") {
            REQUIRE( tree->getNGenNodes() == 0 );
            REQUIRE( tree->getNAllocGenNodes() == 0 );
        }
        WHEN("a non-existing node is fetched") {
            MWNode<3> &node = tree->getNode(r, depth);
            THEN("there will be allocated GenNodes") {
                REQUIRE( tree->getNGenNodes() > 0 );
                REQUIRE( tree->getNAllocGenNodes() > 0 );
            }
            AND_WHEN("the GenNodes are cleared") {
                tree->clearGenerated();
                THEN("there will be un-allocated GenNodes") {
                    REQUIRE( tree->getNGenNodes() > 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
            AND_WHEN("the GenNodes are deleted") {
                tree->deleteGenerated();
                THEN("there will be no GenNodes") {
                    REQUIRE( tree->getNGenNodes() == 0 );
                    REQUIRE( tree->getNAllocGenNodes() == 0 );
                }
            }
        }
        finalize(&tree);
    }
}

} // namespace

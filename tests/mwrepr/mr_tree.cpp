#include "catch.hpp"

#include "factory_functions.h"

namespace mr_tree {

template<int D> void testConstructors();
template<int D> void testNodeFetchers();

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

TEST_CASE("Fetching nodes", "[fetch_nodes], [mr_tree], [trees]") {
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

    const FunctionTree<D> *const_tree = const_cast<const FunctionTree<D> *>(tree);

    // Fetch by NodeIndex
    SECTION("Find node by NodeIndex: existing node") {
        MRNode<D> *node = tree->findNode(idx_0);
        REQUIRE( node != 0 );
        REQUIRE( node->getNodeIndex() == idx_0 );
    }
    SECTION("Find const node by NodeIndex: existing node") {
        const MRNode<D> *node = const_tree->findNode(idx_0);
        REQUIRE( node != 0 );
        REQUIRE( node->getNodeIndex() == idx_0 );
    }
    SECTION("Find node by NodeIndex: non-existing node") {
        MRNode<D> *node = tree->findNode(idx_1);
        REQUIRE( node == 0 );
    }
    SECTION("Find const node by NodeIndex: non-existing node") {
        const MRNode<D> *node = const_tree->findNode(idx_2);
        REQUIRE( node == 0 );
    }
    SECTION("Get node by NodeIndex: existing node") {
        MRNode<D> &node = tree->getNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
    }
    SECTION("Get node by NodeIndex: non-existing node") {
        MRNode<D> &node = tree->getNode(idx_2);
        REQUIRE( node.getNodeIndex() == idx_2 );
    }
    SECTION("Get node or end node by NodeIndex: existing node") {
        MRNode<D> &node = tree->getNodeOrEndNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
    }
    SECTION("Get const node or end node by NodeIndex: existing node") {
        const MRNode<D> &node = const_tree->getNodeOrEndNode(idx_0);
        REQUIRE( node.getNodeIndex() == idx_0 );
    }
    SECTION("Get node or end node by NodeIndex: non-existing node") {
        MRNode<D> &node = tree->getNodeOrEndNode(idx_2);
        REQUIRE( node.getNodeIndex() != idx_0 );
        REQUIRE( node.isEndNode() );
    }
    SECTION("Get const node or end node by NodeIndex: non-existing node") {
        const MRNode<D> &node = const_tree->getNodeOrEndNode(idx_1);
        REQUIRE( node.getNodeIndex() != idx_0 );
        REQUIRE( node.isEndNode() );
    }

    // Fetch by coordinate
    SECTION("Get node by coord: existing node") {
        int depth = 0;
        MRNode<D> &node = tree->getNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node by coord: non-existing node") {
        int depth = 3;
        MRNode<D> &node = tree->getNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node or end node by coord: existing node") {
        int depth = 0;
        MRNode<D> &node = tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get const node or end node by coord: existing node") {
        int depth = 0;
        const MRNode<D> &node = const_tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.getDepth() == depth );
    }
    SECTION("Get node or end node by coord: non-existing node") {
        int depth = 3;
        MRNode<D> &node = tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isEndNode() );
        REQUIRE( node.getDepth() != depth );
    }
    SECTION("Get const node or end node by coord: non-existing node") {
        int depth = 3;
        const MRNode<D> &node = const_tree->getNodeOrEndNode(r, depth);
        REQUIRE( node.hasCoord(r) );
        REQUIRE( node.isEndNode() );
        REQUIRE( node.getDepth() != depth );
    }

    finalize(&tree);
}

} // namespace

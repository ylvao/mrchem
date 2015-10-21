#include "mwtest.h"
#include "BoundingBox.h"
#include "MockTree.h"
#include "MockNode.h"
#include "MathUtils.h"

using namespace std;

class MRTreeTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
        SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
        K = 3;
        tree_1D = 0;
        tree_2D = 0;
        tree_3D = 0;
        initTree<1>(&tree_1D);
        initTree<2>(&tree_2D);
        initTree<3>(&tree_3D);
    }
    virtual void TearDown() {
        deleteTree<1>(&tree_1D);
        deleteTree<2>(&tree_2D);
        deleteTree<3>(&tree_3D);
        ASSERT_TRUE(tree_1D == 0);
        ASSERT_TRUE(tree_2D == 0);
        ASSERT_TRUE(tree_3D == 0);
    }

    template<int D> void initTree(MockTree<D> **tree);
    template<int D> void deleteTree(MockTree<D> **tree);

    template<int D> void initWorldBox(BoundingBox<D> **box);
    template<int D> void deleteWorldBox(BoundingBox<D> **box);

    template<int D> void testTreeInit(const MockTree<D> *tree);

    int K;

    MockTree<1> *tree_1D;
    MockTree<2> *tree_2D;
    MockTree<3> *tree_3D;
};

template<int D>
void MRTreeTest::initWorldBox(BoundingBox<D> **box) {
    int nBoxes[D];
    double origin[D];
    int n = -1;
    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = -(d+1);
        nBoxes[d] = 2*(d+1);
        origin[d] = (double) -d;
    }
    NodeIndex<D> root(n, l);
    *box = new BoundingBox<D>(root, nBoxes, origin);
}

template<int D>
void MRTreeTest::deleteWorldBox(BoundingBox<D> **box) {
    ASSERT_TRUE(box != 0);
    ASSERT_TRUE(*box != 0);
    delete *box;
    *box = 0;
}

template<int D>
void MRTreeTest::initTree(MockTree<D> **tree) {
    ASSERT_TRUE(tree != 0);
    ASSERT_TRUE(*tree == 0);

    BoundingBox<D> *world = 0;
    initWorldBox<D>(&world);
    *tree = new MockTree<D>(K, world);

    deleteWorldBox<D>(&world);
    ASSERT_TRUE(world == 0);
}

template<int D>
void MRTreeTest::deleteTree(MockTree<D> **tree) {
    ASSERT_TRUE(tree != 0);
    ASSERT_TRUE(*tree != 0);
    delete *tree;
    *tree = 0;
}

template<int D>
void MRTreeTest::testTreeInit(const MockTree<D> *tree) {
    ASSERT_TRUE(tree != 0);

    const int k = 3;
    const int kp1 = k+1;
    const int kp1_d = MathUtils::ipow(kp1, D);
    const int tDim = MathUtils::ipow(2, D);

    ASSERT_EQ(k, tree->getOrder());
    ASSERT_EQ(kp1, tree->getKp1());
    ASSERT_EQ(kp1_d, tree->getKp1_d());
    ASSERT_EQ(D, tree->getDim());
    ASSERT_EQ(tDim, tree->getTDim());

    const int rootScale = -1;
    ASSERT_EQ(rootScale, tree->getRootScale());
    ASSERT_EQ(30, tree->getMaxDepth());
    ASSERT_EQ(28, tree->getMaxScale());

    int totNodes = 1;
    for (int d = 0; d < D; d++) {
        const int l = -(d+1);
        const int nodes = 2*(d+1);
        totNodes *= nodes;

        const double unit = pow(2.0, -rootScale);
        const double origin = (double) -d;
        const double length = nodes*unit;
        const double lower = l*unit;
        const double upper = lower + length;

        const double originError = fabs(origin - tree->getOrigin()[d]);
        const double lowerError = fabs(lower - tree->getLowerBounds()[d]);
        const double upperError = fabs(upper - tree->getUpperBounds()[d]);

        EXPECT_LT(originError, MachineZero);
        EXPECT_LT(lowerError, MachineZero);
        EXPECT_LT(upperError, MachineZero);
    }

    ASSERT_EQ(totNodes, tree->getNNodes(-1));
    ASSERT_EQ(totNodes, tree->getNNodes(0));
    ASSERT_EQ(0, tree->getNNodes(1));
    ASSERT_EQ(totNodes, tree->getNEndNodes());
    ASSERT_EQ(totNodes, tree->getNRootNodes());
}

TEST_F(MRTreeTest, Initialization) {
    testTreeInit<1>(tree_1D);
    testTreeInit<2>(tree_2D);
    testTreeInit<3>(tree_3D);
}

TEST_F(MRTreeTest, CopyConstructor) {
    ASSERT_TRUE(tree_1D != 0);
    MockTree<1> *copy_1D = new MockTree<1>(*tree_1D);
    testTreeInit<1>(copy_1D);
    deleteTree<1>(&copy_1D);
    ASSERT_TRUE(copy_1D == 0);

    ASSERT_TRUE(tree_2D != 0);
    MockTree<2> *copy_2D = new MockTree<2>(*tree_2D);
    testTreeInit<2>(copy_2D);
    deleteTree<2>(&copy_2D);
    ASSERT_TRUE(copy_2D == 0);

    ASSERT_TRUE(tree_3D != 0);
    MockTree<3> *copy_3D = new MockTree<3>(*tree_3D);
    testTreeInit<3>(copy_3D);
    deleteTree<3>(&copy_3D);
    ASSERT_TRUE(copy_3D == 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}





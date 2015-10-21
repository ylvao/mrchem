#include "mwtest.h"
#include "MockNode.h"
#include "MockTree.h"

using namespace std;

class MRNodeTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
        SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
        node_1D = 0;
        node_2D = 0;
        node_3D = 0;
        initNode<1>(&node_1D);
        initNode<2>(&node_2D);
        initNode<3>(&node_3D);
    }
    virtual void TearDown() {
        deleteNode<1>(&node_1D);
        deleteNode<2>(&node_2D);
        deleteNode<3>(&node_3D);
        ASSERT_TRUE(node_1D == 0);
        ASSERT_TRUE(node_2D == 0);
        ASSERT_TRUE(node_3D == 0);
    }

    template<int D> void initNode(MockNode<D> **node);
    template<int D> void deleteNode(MockNode<D> **node);

    template<int D> void testNodeInit(const MockNode<D> *node);

    MockNode<1> *node_1D;
    MockNode<2> *node_2D;
    MockNode<3> *node_3D;
};

template<int D>
void MRNodeTest::initNode(MockNode<D> **node) {
    ASSERT_TRUE(node != 0);
    ASSERT_TRUE(*node == 0);
}

template<int D>
void MRNodeTest::deleteNode(MockNode<D> **node) {
    ASSERT_TRUE(node != 0);
    ASSERT_TRUE(*node != 0);
}

template<int D>
void MRNodeTest::testNodeInit(const MockNode<D> *node) {
    ASSERT_TRUE(node != 0);
}

TEST_F(MRNodeTest, Initialization) {
//    testNodeInit<1>(node_1D);
//    testNodeInit<2>(node_2D);
//    testNodeInit<3>(node_3D);
}

TEST_F(MRNodeTest, DefaultConstructor) {
}

TEST_F(MRNodeTest, TreeConstructor) {
}

TEST_F(MRNodeTest, ParentConstructor) {
}

TEST_F(MRNodeTest, CopyConstructorTree) {
//    ASSERT_TRUE(node_1D != 0);
//    MockNode<1> *copy_1D = new MockNode<1>(*node_1D);
//    testNodeInit<1>(copy_1D);
//    deleteNode<1>(&copy_1D);
//    ASSERT_TRUE(copy_1D == 0);

//    ASSERT_TRUE(node_2D != 0);
//    MockNode<2> *copy_2D = new MockNode<2>(*node_2D);
//    testNodeInit<2>(copy_2D);
//    deleteNode<2>(&copy_2D);
//    ASSERT_TRUE(copy_2D == 0);

//    ASSERT_TRUE(node_3D != 0);
//    MockNode<3> *copy_3D = new MockNode<3>(*node_3D);
//    testNodeInit<3>(copy_3D);
//    deleteNode<3>(&copy_3D);
//    ASSERT_TRUE(copy_3D == 0);
}

TEST_F(MRNodeTest, CopyConstructorParent) {
//    ASSERT_TRUE(node_1D != 0);
//    MockNode<1> *copy_1D = new MockNode<1>(*node_1D);
//    testNodeInit<1>(copy_1D);
//    deleteNode<1>(&copy_1D);
//    ASSERT_TRUE(copy_1D == 0);

//    ASSERT_TRUE(node_2D != 0);
//    MockNode<2> *copy_2D = new MockNode<2>(*node_2D);
//    testNodeInit<2>(copy_2D);
//    deleteNode<2>(&copy_2D);
//    ASSERT_TRUE(copy_2D == 0);

//    ASSERT_TRUE(node_3D != 0);
//    MockNode<3> *copy_3D = new MockNode<3>(*node_3D);
//    testNodeInit<3>(copy_3D);
//    deleteNode<3>(&copy_3D);
//    ASSERT_TRUE(copy_3D == 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}




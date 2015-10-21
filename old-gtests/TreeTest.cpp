#include "mwtest.h"
#include "NodeBox.h"
#include "MRGrid.h"
#include "GridNode.h"
#include "NodeIndex.h"
#include "GridGenerator.h"
#include "LebesgueIterator.h"

using namespace std;

class TreeTest: public ::testing::Test {
public:
    static void SetUpTestCase() {
	SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }
};

TEST_F(TreeTest, Initialization) {
    int rootScale = -2;
    int nBoxes[1] = {3};
    double origin[1] = {6.0};
    NodeBox<1> box(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<1> grid(order, &box);
    const double *lb = grid.getLowerBounds();
    const double *ub = grid.getUpperBounds();
    const double *l = grid.getRootBox().getBoxLength();

    EXPECT_EQ(grid.getNNodes(), 3);
    EXPECT_EQ(grid.getNNodes(0), 3);
    EXPECT_EQ(grid.getNNodes(1), 0);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
}

TEST_F(TreeTest, GetNodeByIndex) {
    int rootScale = -2;
    int nBoxes[1] = {3};
    double origin[1] = {6.0};
    NodeBox<1> box(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<1> grid(order, &box);
    const MRGrid<1> &constGrid = const_cast<MRGrid<1> &>(grid);

    int n = 0;
    int l[1] = {0};
    NodeIndex<1> idx(n, l);
    const MRNode<1> *constNode = constGrid.findNode(idx);
    int pt = -1;
    if (constNode == 0) pt = 0;
    EXPECT_EQ(pt, 0);
    EXPECT_EQ(grid.getNNodes(), 3);
    EXPECT_EQ(grid.getNEndNodes(), 3);

    const MRNode<1> *node_1 = grid.findNode(idx);
    pt = -1;
    if (node_1 == 0) pt = 0;
    EXPECT_EQ(pt, 0);
    EXPECT_EQ(grid.getNNodes(), 3);
    EXPECT_EQ(grid.getNEndNodes(), 3);
    EXPECT_EQ(grid.getNGenNodes(), 0);
    EXPECT_EQ(grid.getNAllocGenNodes(), 0);

    const MRNode<1> &node_2 = grid.getNode(idx);
    EXPECT_EQ(node_2.getScale(), n);
    EXPECT_EQ(node_2.getTranslation()[0], l[0]);
    EXPECT_EQ(grid.getNNodes(), 7);
    EXPECT_EQ(grid.getNEndNodes(), 3);
    EXPECT_EQ(grid.getNGenNodes(), 4);
    EXPECT_EQ(grid.getNAllocGenNodes(), 4);

    grid.purgeGenerated();
    EXPECT_EQ(grid.getNNodes(), 3);
    EXPECT_EQ(grid.getNEndNodes(), 3);
    EXPECT_EQ(grid.getNGenNodes(), 0);
    EXPECT_EQ(grid.getNAllocGenNodes(), 0);

    idx.setScale(rootScale);
    const MRNode<1> *node_3 = grid.findNode(idx);
    pt = -1;
    if (node_3 == 0) pt = 0;
    EXPECT_NE(pt, 0);
    EXPECT_EQ(node_3->getScale(), rootScale);
    EXPECT_EQ(node_3->getTranslation()[0], l[0]);
    EXPECT_EQ(grid.getNNodes(), 3);
    EXPECT_EQ(grid.getNEndNodes(), 3);
    EXPECT_EQ(grid.getNGenNodes(), 0);
    EXPECT_EQ(grid.getNAllocGenNodes(), 0);
}

TEST_F(TreeTest, GetNodeByCoord) {
    double r0[3] = {0.0,0.0,0.0};
    double r1[3] = {pi/10.0,pi/10.0,pi/10.0};

    int rootScale = 0;
    int nBoxes[3] = {2,2,2};
    double origin[3] = {1.0,1.0,1.0};
    NodeBox<3> box(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<3> grid(order, &box);

    GridGenerator<3> generator;
    generator.setUniformScale(3);
    generator.generateGrid(grid);

    println(0, endl);
    println(0, "const find");

    const MRGrid<3> &constGrid = const_cast<MRGrid<3> &>(grid);
    const MRNode<3> *constNode = constGrid.findNode(r0);
    EXPECT_EQ(constNode->getScale(), 3);
    EXPECT_EQ(constNode->getTranslation()[0], 0);
    EXPECT_EQ(constNode->getTranslation()[1], 0);
    EXPECT_EQ(constNode->getTranslation()[2], 0);

    println(0, endl);
    println(0, "non-const find");

    MRNode<3> *node_1 = grid.findNode(r0);
    EXPECT_EQ(node_1->getScale(), 3);
    EXPECT_EQ(node_1->getTranslation()[0], 0);
    EXPECT_EQ(node_1->getTranslation()[1], 0);
    EXPECT_EQ(node_1->getTranslation()[2], 0);

    MRNode<3> *node_2 = grid.findNode(r0, 1);
    EXPECT_EQ(node_2->getScale(), 1);
    EXPECT_EQ(node_2->getTranslation()[0], 0);
    EXPECT_EQ(node_2->getTranslation()[1], 0);
    EXPECT_EQ(node_2->getTranslation()[2], 0);

    println(0, endl);

    MRNode<3> *node_3 = grid.findNode(r1, 4);
    EXPECT_EQ(node_2->getScale(), 3);
    EXPECT_EQ(node_2->getTranslation()[0], 0);
    EXPECT_EQ(node_2->getTranslation()[1], 0);
    EXPECT_EQ(node_2->getTranslation()[2], 0);
}

/*
TEST_F(TreeTest, Distribution) {
    int rootScale = 0;
    int nBoxes[3] = {1,1,2};
    double origin[3] = {0.0, 0.0, 1.0};
    NodeBox<3> box(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<3> grid(order, &box);

    println(0, grid);

    GridGenerator<3> generator;
    generator.setUniformScale(1);

    generator.generateGrid(grid);

    println(0, grid);
    grid.printNodeRankCount();

    LebesgueIterator<3> it_1(&grid);
    while(it_1.next()) {
	MRNode<3> &node = it_1.getNode();
	println(0, node);
    }
    println(0, endl);

    grid.distributeEndNodes();
    grid.printNodeRankCount();

    LebesgueIterator<3> it_2(&grid);
    while(it_2.next()) {
	MRNode<3> &node = it_2.getNode();
	println(0, node);
    }
    println(0, endl);


    return;
}
*/
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::AddGlobalTestEnvironment(new Environment);
    mpi::environment env(argc, argv);
    //MREnv::initializeMRCPP(argc, argv, "FuncTreeTest");
    return RUN_ALL_TESTS();
}



#include "mwtest.h"
#include "BoundingBox.h"

using namespace std;

class BoundingBoxTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
        SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
        box_1D = 0;
        box_2D = 0;
        box_3D = 0;
        initBox<1>(&box_1D);
        initBox<2>(&box_2D);
        initBox<3>(&box_3D);
    }
    virtual void TearDown() {
        deleteBox<1>(&box_1D);
        deleteBox<2>(&box_2D);
        deleteBox<3>(&box_3D);
        ASSERT_TRUE(box_1D == 0);
        ASSERT_TRUE(box_2D == 0);
        ASSERT_TRUE(box_3D == 0);
    }

    template<int D> void initBox(BoundingBox<D> **box);
    template<int D> void deleteBox(BoundingBox<D> **box);

    template<int D> void testBoxInit(const BoundingBox<D> *box);
    template<int D> void testNodeIndexByCoord(const BoundingBox<D> *box);
    template<int D> void testNodeIndexByBox(const BoundingBox<D> *box);
    template<int D> void testBoxIndexByCoord(const BoundingBox<D> *box);
    template<int D> void testBoxIndexByNode(const BoundingBox<D> *box);

    BoundingBox<1> *box_1D;
    BoundingBox<2> *box_2D;
    BoundingBox<3> *box_3D;
};

template<int D>
void BoundingBoxTest::initBox(BoundingBox<D> **box) {
    ASSERT_TRUE(box != 0);
    ASSERT_TRUE(*box == 0);
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
void BoundingBoxTest::deleteBox(BoundingBox<D> **box) {
    ASSERT_TRUE(box != 0);
    ASSERT_TRUE(*box != 0);
    delete *box;
    *box = 0;
}

template<int D>
void BoundingBoxTest::testBoxInit(const BoundingBox<D> *box) {
    ASSERT_TRUE(box != 0);

    const int scale = -1;
    EXPECT_EQ(scale, box->getRootScale());

    const double unit = pow(2.0, -scale);
    const double unitError = fabs(unit - box->getUnitLength());
    EXPECT_LT(unitError, MachineZero);

    int totBoxes = 1;
    for (int d = 0; d < D; d++) {
        const int l = -(d+1);
        const int boxes = 2*(d+1);

        EXPECT_EQ(l, box->getCornerIndex().getTranslation()[d]);
        EXPECT_EQ(boxes, box->getNBoxes(d));

        const double origin = (double) -d;
        const double length = boxes*unit;
        const double lower = l*unit;
        const double upper = lower + length;

        const double originError = fabs(origin - box->getOrigin()[d]);
        const double lengthError = fabs(length - box->getBoxLength()[d]);
        const double lowerError = fabs(lower - box->getLowerBounds()[d]);
        const double upperError = fabs(upper - box->getUpperBounds()[d]);

        EXPECT_LT(originError, MachineZero);
        EXPECT_LT(lengthError, MachineZero);
        EXPECT_LT(lowerError, MachineZero);
        EXPECT_LT(upperError, MachineZero);

        totBoxes *= boxes;
    }
    EXPECT_EQ(totBoxes, box->getNBoxes());
}

template<int D>
void BoundingBoxTest::testNodeIndexByBox(const BoundingBox<D> *box) {
    ASSERT_TRUE(box != 0);

    const int n = -1;
    const int L1[3] = {0,0,0};
    const int L2[3] = {-1,-2,2};
    const NodeIndex<D> ref1(n,L1);
    const NodeIndex<D> ref2(n,L2);

    if (D == 1) {
        EXPECT_TRUE(ref1 == box->getNodeIndex(1));
        EXPECT_TRUE(ref2 == box->getNodeIndex(0));
    }
    if (D == 2) {
        EXPECT_TRUE(ref1 == box->getNodeIndex(5));
        EXPECT_TRUE(ref2 == box->getNodeIndex(0));
    }
    if (D == 3) {
        EXPECT_TRUE(ref1 == box->getNodeIndex(29));
        EXPECT_TRUE(ref2 == box->getNodeIndex(40));
    }
}

template<int D>
void BoundingBoxTest::testNodeIndexByCoord(const BoundingBox<D> *box) {
    ASSERT_TRUE(box != 0);

    const int n = -1;
    const double a = -pi/10.0;
    const double r1[3] = {0.0,0.0,0.0};
    const double r2[3] = {a + 0.0, a + 1.0, a + 2.0};

    const int L1[3] = {0,0,0};
    const NodeIndex<D> ref1(n, L1);
    EXPECT_TRUE(ref1 == box->getNodeIndex(r1));

    const int L2[3] = {-1,0,0};
    const NodeIndex<D> ref2(n, L2);
    EXPECT_TRUE(ref2 == box->getNodeIndex(r2));
}

template<int D>
void BoundingBoxTest::testBoxIndexByNode(const BoundingBox<D> *box) {
    ASSERT_TRUE(box != 0);
    const int n1 = -1;
    const int n2 = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
        l[d] = d*(1 - d) - 1;
    }

    const NodeIndex<D> idx_1(n1, l);
    const NodeIndex<D> idx_2(n2, l);

    if (D == 1) {
        EXPECT_EQ(0, box->getBoxIndex(idx_1));
        EXPECT_EQ(0, box->getBoxIndex(idx_2));
    }
    if (D == 2) {
        EXPECT_EQ(2, box->getBoxIndex(idx_1));
        EXPECT_EQ(2, box->getBoxIndex(idx_2));
    }
    if (D == 3) {
        EXPECT_EQ(2, box->getBoxIndex(idx_1));
        EXPECT_EQ(18, box->getBoxIndex(idx_2));
    }
}

template<int D>
void BoundingBoxTest::testBoxIndexByCoord(const BoundingBox<D> *box) {
    ASSERT_TRUE(box != 0);
    double r1[D];
    double r2[D];
    for (int d = 0; d < D; d++) {
        r1[d] = 0.0;
        r2[d] = -pi/10.0 + d;
    }
    if (D == 1) {
        EXPECT_EQ(1, box->getBoxIndex(r1));
        EXPECT_EQ(0, box->getBoxIndex(r2));
    }
    if (D == 2) {
        EXPECT_EQ(5, box->getBoxIndex(r1));
        EXPECT_EQ(4, box->getBoxIndex(r2));
    }
    if (D == 3) {
        EXPECT_EQ(29, box->getBoxIndex(r1));
        EXPECT_EQ(28, box->getBoxIndex(r2));
    }
}

TEST_F(BoundingBoxTest, Initialization) {
    testBoxInit<1>(box_1D);
    testBoxInit<2>(box_2D);
    testBoxInit<3>(box_3D);
}


TEST_F(BoundingBoxTest, CopyConstructor) {
    ASSERT_TRUE(box_1D != 0);
    BoundingBox<1> *copy_1D = new BoundingBox<1>(*box_1D);
    testBoxInit<1>(copy_1D);
    deleteBox<1>(&copy_1D);
    ASSERT_TRUE(copy_1D == 0);

    ASSERT_TRUE(box_2D != 0);
    BoundingBox<2> *copy_2D = new BoundingBox<2>(*box_2D);
    testBoxInit<2>(copy_2D);
    deleteBox<2>(&copy_2D);
    ASSERT_TRUE(copy_2D == 0);

    ASSERT_TRUE(box_3D != 0);
    BoundingBox<3> *copy_3D = new BoundingBox<3>(*box_3D);
    testBoxInit<3>(copy_3D);
    deleteBox<3>(&copy_3D);
    ASSERT_TRUE(copy_3D == 0);
}

TEST_F(BoundingBoxTest, GetNodeIndexByCoord) {
    testNodeIndexByCoord<1>(box_1D);
    testNodeIndexByCoord<2>(box_2D);
    testNodeIndexByCoord<3>(box_3D);
}

TEST_F(BoundingBoxTest, GetNodeIndexByBox) {
    testNodeIndexByBox<1>(box_1D);
    testNodeIndexByBox<2>(box_2D);
    testNodeIndexByBox<3>(box_3D);
}

TEST_F(BoundingBoxTest, GetBoxIndexByCoord) {
    testBoxIndexByCoord<1>(box_1D);
    testBoxIndexByCoord<2>(box_2D);
    testBoxIndexByCoord<3>(box_3D);
}

TEST_F(BoundingBoxTest, GetBoxIndexByNode) {
    testBoxIndexByNode<1>(box_1D);
    testBoxIndexByNode<2>(box_2D);
    testBoxIndexByNode<3>(box_3D);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}



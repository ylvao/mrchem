#include "mwtest.h"
#include "NodeBox.h"

using namespace std;

class NodeBoxTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
        SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
        this->box_1D = 0;
        this->box_2D = 0;
        this->box_3D = 0;
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

    template<int D> void initBox(NodeBox<D> **box);
    template<int D> void deleteBox(NodeBox<D> **box);

    template<int D> void testBoxInit(const NodeBox<D> *box);

    NodeBox<1> *box_1D;
    NodeBox<2> *box_2D;
    NodeBox<3> *box_3D;
};

template<int D>
void NodeBoxTest::initBox(NodeBox<D> **box) {
    ASSERT_TRUE(box != 0);
    ASSERT_TRUE(*box == 0);
    int n = -1;
    int l[D];
    int nBoxes[D];
    double origin[D];
    for (int d = 0; d < D; d++) {
        l[d] = -(d+1);
        nBoxes[d] = 2*(d+1);
        origin[d] = (double) -d;
    }
    NodeIndex<D> root(n, l);
    *box = new NodeBox<D>(root, nBoxes, origin);
}

template<int D>
void NodeBoxTest::deleteBox(NodeBox<D> **box) {
    ASSERT_TRUE(box != 0);
    ASSERT_TRUE(*box != 0);
    delete *box;
    *box = 0;
}

template<int D>
void NodeBoxTest::testBoxInit(const NodeBox<D> *box) {
    ASSERT_TRUE(box != 0);
}

TEST_F(NodeBoxTest, NodeBoxInit) {
    testBoxInit<1>(box_1D);
    testBoxInit<2>(box_2D);
    testBoxInit<3>(box_3D);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}




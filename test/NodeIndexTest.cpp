#include "mwtest.h"
#include "NodeIndex.h"

using namespace std;

class NodeIndexTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
	SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
	idx_1D = 0;
	idx_2D = 0;
	idx_3D = 0;
	initializeIndex<1>(&idx_1D);
	initializeIndex<2>(&idx_2D);
	initializeIndex<3>(&idx_3D);
    }

    virtual void TearDown() {
	deleteIndex<1>(&idx_1D);
	deleteIndex<2>(&idx_2D);
	deleteIndex<3>(&idx_3D);
	ASSERT_TRUE(idx_1D == 0);
	ASSERT_TRUE(idx_2D == 0);
	ASSERT_TRUE(idx_3D == 0);
    }

    template<int D> void initializeIndex(NodeIndex<D> **idx);
    template<int D> void deleteIndex(NodeIndex<D> **idx);

    template<int D> void initializeScale(NodeIndex<D> **idx, int n);
    template<int D> void initializeRank(NodeIndex<D> **idx, int r);
    template<int D> void initializeTranslation(NodeIndex<D> **idx, const int *l);

    template<int D> void testInitialization(const NodeIndex<D> *idx);

    NodeIndex<1> *idx_1D;
    NodeIndex<2> *idx_2D;
    NodeIndex<3> *idx_3D;
};

template<int D> 
void NodeIndexTest::initializeIndex(NodeIndex<D> **idx) {
    ASSERT_TRUE(idx != 0);
    ASSERT_TRUE(*idx == 0);
    int n = 1;
    int r = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
	l[d] = d-1;
    }
    *idx = new NodeIndex<D>(n, l, r);
}

template<int D> 
void NodeIndexTest::deleteIndex(NodeIndex<D> **idx) {
    ASSERT_TRUE(idx != 0);
    ASSERT_TRUE(*idx != 0);
    delete *idx;
    *idx = 0;
}

template<int D> 
void NodeIndexTest::initializeScale(NodeIndex<D> **idx, int n) {
    ASSERT_TRUE(idx != 0);
    ASSERT_TRUE(*idx == 0);
    int r = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
	l[d] = d-1;
    }
    *idx = new NodeIndex<D>(n, l, r);
}

template<int D> 
void NodeIndexTest::initializeRank(NodeIndex<D> **idx, int r) {
    ASSERT_TRUE(idx != 0);
    ASSERT_TRUE(*idx == 0);
    int n = 1;
    int l[D];
    for (int d = 0; d < D; d++) {
	l[d] = d-1;
    }
    *idx = new NodeIndex<D>(n, l, r);
}

template<int D> 
void NodeIndexTest::initializeTranslation(NodeIndex<D> **idx, const int *l) {
    ASSERT_TRUE(idx != 0);
    ASSERT_TRUE(*idx == 0);
    int n = 1;
    int r = 1;
    *idx = new NodeIndex<D>(n, l, r);
}

template<int D>
void NodeIndexTest::testInitialization(const NodeIndex<D> *idx) {
    ASSERT_TRUE(idx != 0);

    const int scale = 1;
    EXPECT_EQ(scale, idx->getScale());

    const int rank = 1;
    EXPECT_EQ(1, idx->getRankId());

    for (int d = 0; d < D; d++) {
        const int l = d-1;
	EXPECT_EQ(l, idx->getTranslation(d));
	EXPECT_EQ(l, idx->getTranslation()[d]);
    }
}

TEST_F(NodeIndexTest, constructor) {
    testInitialization<1>(idx_1D);
    testInitialization<2>(idx_2D);
    testInitialization<3>(idx_3D);
}

TEST_F(NodeIndexTest, copyContructor_1D) {
    ASSERT_TRUE(idx_1D != 0);
    NodeIndex<1> *copy_1D = new NodeIndex<1>(*idx_1D);
    testInitialization(copy_1D);
    deleteIndex(&copy_1D);
    ASSERT_TRUE(copy_1D == 0);
}

TEST_F(NodeIndexTest, copyContructor_2D) {
    ASSERT_TRUE(idx_2D != 0);
    NodeIndex<2> *copy_2D = new NodeIndex<2>(*idx_2D);
    testInitialization(copy_2D);
    deleteIndex(&copy_2D);
    ASSERT_TRUE(copy_2D == 0);
}

TEST_F(NodeIndexTest, copyContructor_3D) {
    ASSERT_TRUE(idx_3D != 0);
    NodeIndex<3> *copy_3D = new NodeIndex<3>(*idx_3D);
    testInitialization(copy_3D);
    deleteIndex(&copy_3D);
    ASSERT_TRUE(copy_3D == 0);
}

TEST_F(NodeIndexTest, assignmentOperator_1D) {
    ASSERT_TRUE(idx_1D != 0);
    NodeIndex<1> *copy_1D = new NodeIndex<1>();
    *copy_1D = *idx_1D;
    testInitialization(copy_1D);
    deleteIndex(&copy_1D);
    ASSERT_TRUE(copy_1D == 0);
}

TEST_F(NodeIndexTest, assignmentOperator_2D) {
    ASSERT_TRUE(idx_2D != 0);
    NodeIndex<2> *copy_2D = new NodeIndex<2>();
    *copy_2D = *idx_2D;
    testInitialization(copy_2D);
    deleteIndex(&copy_2D);
    ASSERT_TRUE(copy_2D == 0);
}

TEST_F(NodeIndexTest, assignmentOperator_3D) {
    ASSERT_TRUE(idx_3D != 0);
    NodeIndex<3> *copy_3D = new NodeIndex<3>();
    *copy_3D = *idx_3D;
    testInitialization(copy_3D);
    deleteIndex(&copy_3D);
    ASSERT_TRUE(copy_3D == 0);
}

TEST_F(NodeIndexTest, compareScale_1D) {
    NodeIndex<1> *idx = 0;
    for (int n = -2; n < 2; n++) {
	ASSERT_TRUE(idx == 0);
        initializeScale<1>(&idx, n);
	if (n == 1) {
	    EXPECT_TRUE(*idx == *idx_1D);
	    EXPECT_FALSE(*idx != *idx_1D);
	} else {
	    EXPECT_FALSE(*idx == *idx_1D);
	    EXPECT_TRUE(*idx != *idx_1D);
	}
	deleteIndex<1>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareScale_2D) {
    NodeIndex<2> *idx = 0;
    for (int n = -2; n < 2; n++) {
	ASSERT_TRUE(idx == 0);
        initializeScale<2>(&idx, n);
	if (n == 1) {
	    EXPECT_TRUE(*idx == *idx_2D);
	    EXPECT_FALSE(*idx != *idx_2D);
	} else {
	    EXPECT_FALSE(*idx == *idx_2D);
	    EXPECT_TRUE(*idx != *idx_2D);
	}
	deleteIndex<2>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareScale_3D) {
    NodeIndex<3> *idx = 0;
    for (int n = -2; n < 2; n++) {
	ASSERT_TRUE(idx == 0);
        initializeScale<3>(&idx, n);
	if (n == 1) {
	    EXPECT_TRUE(*idx == *idx_3D);
	    EXPECT_FALSE(*idx != *idx_3D);
	} else {
	    EXPECT_FALSE(*idx == *idx_3D);
	    EXPECT_TRUE(*idx != *idx_3D);
	}
	deleteIndex<3>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareRank_1D) {
    NodeIndex<1> *idx = 0;
    for (int r = -2; r < 2; r++) {
	ASSERT_TRUE(idx == 0);
        initializeRank<1>(&idx, r);
	EXPECT_TRUE(*idx == *idx_1D);
	EXPECT_FALSE(*idx != *idx_1D);
	deleteIndex<1>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareRank_2D) {
    NodeIndex<2> *idx = 0;
    for (int r = -2; r < 2; r++) {
	ASSERT_TRUE(idx == 0);
        initializeRank<2>(&idx, r);
	EXPECT_TRUE(*idx == *idx_2D);
	EXPECT_FALSE(*idx != *idx_2D);
	deleteIndex<2>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareRank_3D) {
    NodeIndex<3> *idx = 0;
    for (int r = -2; r < 2; r++) {
	ASSERT_TRUE(idx == 0);
        initializeRank<3>(&idx, r);
	EXPECT_TRUE(*idx == *idx_3D);
	EXPECT_FALSE(*idx != *idx_3D);
	deleteIndex<3>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareTranslation_1D) {
    int l[1];
    l[0] = idx_2D->getTranslation(0);
    NodeIndex<1> *idx = 0;
    for (int x = -2; x < 2; x++) {
	ASSERT_TRUE(idx == 0);
	l[0] = x;
        initializeTranslation<1>(&idx, l);
	if (x == -1) {
	    EXPECT_TRUE(*idx == *idx_1D);
	    EXPECT_FALSE(*idx != *idx_1D);
	} else {
	    EXPECT_FALSE(*idx == *idx_1D);
	    EXPECT_TRUE(*idx != *idx_1D);
	}
	deleteIndex<1>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareTranslation_2D) {
    int l[2];
    l[0] = idx_2D->getTranslation(0);
    l[1] = idx_2D->getTranslation(1);
    NodeIndex<2> *idx = 0;
    for (int y = -2; y < 2; y++) {
	ASSERT_TRUE(idx == 0);
	l[1] = y;
        initializeTranslation<2>(&idx, l);
	if (y == 0) {
	    EXPECT_TRUE(*idx == *idx_2D);
	    EXPECT_FALSE(*idx != *idx_2D);
	} else {
	    EXPECT_FALSE(*idx == *idx_2D);
	    EXPECT_TRUE(*idx != *idx_2D);
	}
	deleteIndex<2>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

TEST_F(NodeIndexTest, compareTranslation_3D) {
    int l[3];
    l[0] = idx_3D->getTranslation(0);
    l[1] = idx_3D->getTranslation(1);
    l[2] = idx_3D->getTranslation(2);
    NodeIndex<3> *idx = 0;
    for (int z = -2; z < 2; z++) {
	ASSERT_TRUE(idx == 0);
	l[2] = z;
        initializeTranslation<3>(&idx, l);
	if (z == 1) {
	    EXPECT_TRUE(*idx == *idx_3D);
	    EXPECT_FALSE(*idx != *idx_3D);
	} else {
	    EXPECT_FALSE(*idx == *idx_3D);
	    EXPECT_TRUE(*idx != *idx_3D);
	}
	deleteIndex<3>(&idx);
    }
    ASSERT_TRUE(idx == 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}



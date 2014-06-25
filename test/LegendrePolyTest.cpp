#include "mwtest.h"
#include "LegendrePoly.h"

using namespace std;

class LegendrePolyTest: public ::testing::Test {
protected:
    static void SetUpTestCase() {
	SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
	ASSERT_EQ(0, L.size());
	for (int k = 0; k < 10; k++) {
	    LegendrePoly *leg = new LegendrePoly(k, 2.0, 1.0);
	    L.push_back(leg);
	}
	ASSERT_EQ(10, L.size());
    }

    virtual void TearDown() {
	for (int k = 0; k < L.size(); k++) {
	    if (L[k] != 0) {
		delete L[k];
		L[k] = 0;
	    }
	}
    }

    vector<LegendrePoly *> L;
};

TEST_F(LegendrePolyTest, Initialization) {
    for (int k = 0; k < L.size(); k++) {
	EXPECT_EQ(L[k]->size(), k + 1);

	double lb = L[k]->getLowerBound();
	double ub = L[k]->getUpperBound();
	EXPECT_LT(fabs(lb + 0.0), MachineZero); 
	EXPECT_LT(fabs(ub - 1.0), MachineZero); 

	double x = 1.0;
	double val = L[k]->evalf(&x);
	EXPECT_LT(fabs(val - 1.0), MachineZero); 
    }
}

TEST_F(LegendrePolyTest, Orthogonality) {
    ASSERT_GT(L.size(), 0);
    for (int i = 0; i < L.size(); i++) {
	ASSERT_TRUE(L[i] != 0);
	LegendrePoly &L_i = *L[i];
	for (int j = 0; j < L.size(); j++) {
	    ASSERT_TRUE(L[j] != 0);
	    LegendrePoly &L_j = *L[j];
	    const double overlap = L_i.innerProduct(L_j);
	    if (i == j) {
		EXPECT_GT(fabs(overlap), MachineZero);
	    } else {
		EXPECT_LT(fabs(overlap), MachineZero);
	    }
	}
    }
}

    
int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}



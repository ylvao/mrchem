#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"

using namespace std;

class SkeletonTest: public ::testing::Test {
public:
	static void SetUpTestCase() {
		GaussExp<3> f = MRUtils::genGaussExp<3>(1, 0.5, 0.1, 3000.0);
		defFunc = new GaussExp<3>(f);
		defTree = new FunctionTree<3>();
		defTree->projectFunction(*defFunc);
		defTree->setName("defTree");
//		coulomb = new CoulombOperator();
	}
	static void TearDownTestCase() {
		delete defFunc;
		delete defTree;
		delete coulomb;
	}

	virtual void SetUp() {
		SET_PRINT_PRECISION(15);
		cout << scientific << setprecision(15);
		rank = 0;
		mpi::communicator world;
		rank = world.rank();
	}

	virtual void TearDown() {
	}

	int rank;

	static FunctionTree<3> *defTree;
	static GaussExp<3> *defFunc;
	static CoulombOperator *coulomb;
};

TEST_F(SkeletonTest, MyLittleTester) {
	EXPECT_EQ(0.0, 1.0);
}

FunctionTree<3> *SkeletonTest::defTree = 0;
GaussExp<3> *SkeletonTest::defFunc = 0;
CoulombOperator *SkeletonTest::coulomb = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "SkeletonTest");
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}


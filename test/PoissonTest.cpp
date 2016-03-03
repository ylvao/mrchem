#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"

using namespace std;

class PoissonTest: public ::testing::Test {
public:
	static void SetUpTestCase() {
		GaussExp<3> f = MRUtils::genGaussExp<3>(1, 0.5, 0.1, 3000.0);
		defFunc = new GaussExp<3>(f);
		defTree = new FunctionTree<3>();
		defTree->projectFunction(*defFunc);
		defTree->setName("defTree");

		coulomb = new CoulombOperator();
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
		r1[0] = 0.1;
		r1[1] = 0.2;
		r1[2] = 0.3;

		r2[0] = 0.2;
		r2[1] = 0.1;
		r2[2] = 0.0;
	}

	virtual void TearDown() {
	}

	int rank;

	double r1[3];
	double r2[3];

	static FunctionTree<3> *defTree;
	static GaussExp<3> *defFunc;
	static CoulombOperator *coulomb;
};

TEST_F(PoissonTest, ApplyOperator) {
	defTree->setAutoClean(true);

	FunctionTree<3> gTree;
	gTree.setName("gTree");
	coulomb->setUseGenNodes(false);
	gTree = coulomb->apply(*defTree);

	println(0, *defTree);
	println(0, gTree);
	double gint = gTree.integrate();
	println(0, "gTree integral: " << gint);
}

FunctionTree<3> *PoissonTest::defTree = 0;
CoulombOperator *PoissonTest::coulomb = 0;
GaussExp<3> *PoissonTest::defFunc = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "PoissonTest");
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}

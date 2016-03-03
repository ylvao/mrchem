#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "XCFun.h"
#include "DiracExchangeD1.h"

using namespace std;

class XCFunTest: public ::testing::Test {
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

TEST_F(XCFunTest, XCFun) {
	EXPECT_EQ(0, 1);
	return;
	GaussExp<3> f = MRUtils::genGaussExp<3>(1, 1.0, 0.1, 300.0);

	FunctionTree<3> fTree;
	fTree.projectFunction(f);

	XCFun xcFun;
	FunctionTree<3> gTree;
	xcFun.calcXCPotential(fTree, gTree);

	DiracExchangeD1 dirac;
	FunctionTree<3> hTree;
	hTree = dirac(fTree);

	println(0, fTree);
	println(0, gTree);
	println(0, hTree);

	Plot<3> plt;
	double b[3] = {4.0, 4.0, 4.0};
	plt.setUpperBounds(b);
	plt.linePlot(fTree, "fTree");
	plt.linePlot(gTree, "gTree");
	plt.linePlot(hTree, "hTree");

	FunctionTree<3> errorTree;
	errorTree = gTree - hTree;

	double sqError = errorTree.getSquareNorm();
	double gInt = gTree.integrate();
	double hInt = hTree.integrate();

	println(0, "\n\nSquare error: " << setw(21) << sqError/gTree.getSquareNorm());
	println(0, "Old integral:     " << setw(21) << gInt);
	println(0, "XCFun integral:   " << setw(21) << hInt);
	EXPECT_EQ(0.0, 1.0);
}

FunctionTree<3> *XCFunTest::defTree = 0;
GaussExp<3> *XCFunTest::defFunc = 0;
CoulombOperator *XCFunTest::coulomb = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "XCFunTest");
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}



#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "OperatorTree.h"
#include "Derivative.h"
#include "TreeIterator.h"

#define D 3

using namespace std;

class OperatorTest: public ::testing::Test {
public:
	static void SetUpTestCase() {
		GaussExp<D> f = MRUtils::genGaussExp<D>(4, 0.2, 0.2, 300.0);
		defFunc = new GaussExp<D>(f);
		defTree = new FunctionTree<D>();
		defTree->projectFunction(*defFunc);
		defTree->setName("defTree");
	}
	static void TearDownTestCase() {
		delete defFunc;
		delete defTree;
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

	void TestDerivative(Derivative &derv) {
		double error, relError;
		FunctionTree<D> errTree;
		errTree.setRelPrec(errTree.getRelPrec() * 0.01);

		FunctionTree<D> fTree;
		fTree.projectFunction(*defFunc);
		fTree.setName("fTree");
		println(1, fTree);

		for (int d = 0; d < D; d++) {
			FunctionTree<D> dfTree;
			derv.setApplyDir(d);;
			dfTree = derv.apply(fTree);
			dfTree.setName("dfTree");
			println(1, dfTree);

			GaussExp<D> df = defFunc->differentiate(d);

			FunctionTree<D> dFTree;
			dFTree.setRelPrec(dfTree.getRelPrec() * 0.01);
			dFTree.projectFunction(df);
			dFTree.setName("dFTree");
			println(1, dFTree);

			errTree = dfTree - dFTree;
			println(1, errTree);

			error = sqrt(errTree.getSquareNorm());
			relError = error/sqrt(dFTree.getSquareNorm());

			EXPECT_LT(relError, fTree.getRelPrec());
		}
	}

	int rank;

	double r1[3];
	double r2[3];

	static FunctionTree<D> *defTree;
	static GaussExp<D> *defFunc;
};

TEST_F(OperatorTest, TreeDerivative) {
	Derivative derv(0.0, 0.0);
	TestDerivative(derv);
}

FunctionTree<D> *OperatorTest::defTree = 0;
GaussExp<D> *OperatorTest::defFunc = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "OperatorTest");
	FunctionTree<D>::setDefaultPrecision(1.0e-3);
	FunctionTree<D>::setDefaultOrder(7);
	Operator::setDefaultOrder(7);
	OperatorTree::setDefaultOrder(7);
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}





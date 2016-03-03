#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "Plot.h"
#include "TreeIterator.h"

#define D 3

using namespace std;

class TreeArithmetic: public ::testing::Test {
public:
	static void SetUpTestCase() {
		GaussExp<D> f = MRUtils::genGaussExp<D>(1, 0.5, 0.1, 3000.0);
		GaussExp<D> g = MRUtils::genGaussExp<D>(2, 0.45, 0.1, 3000.0);
		defFuncF = new GaussExp<D>(f);
		defFuncG = new GaussExp<D>(g);
		defTreeF = new FunctionTree<D>();
		defTreeG = new FunctionTree<D>();
		defTreeF->projectFunction(*defFuncF);
		defTreeG->projectFunction(*defFuncG);
		defTreeF->setName("defTreeF");
		defTreeG->setName("defTreeG");
	}
	static void TearDownTestCase() {
		delete defFuncF;
		delete defFuncG;
		delete defTreeF;
		delete defTreeG;
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

	void testTrees(FunctionTree<D> &testTree, FunctionTree<D> &refTree);

	int rank;

	double r1[3];
	double r2[3];

	static FunctionTree<D> *defTreeF;
	static FunctionTree<D> *defTreeG;
	static GaussExp<D> *defFuncF;
	static GaussExp<D> *defFuncG;
};

void TreeArithmetic::testTrees(FunctionTree<D> &testTree, FunctionTree<D> &refTree) {
	FunctionTree<D> errTree;
	errTree.add(1.0, testTree, -1.0, refTree);

	double v1 = testTree.evalf(r1);
	double v2 = refTree.evalf(r1);
	double dv = fabs(v1 - v2);

	double errInt = errTree.integrate();
	double relNormErr = sqrt(errTree.getSquareNorm()/testTree.getSquareNorm());

	double testInt = testTree.integrate();
	double refInt = refTree.integrate();
	double diffInt = fabs(testInt - refInt);

	EXPECT_LT(dv, testTree.getRelPrec());
	EXPECT_LT(errInt, testTree.getRelPrec());
	EXPECT_LT(relNormErr, testTree.getRelPrec());
	EXPECT_LT(diffInt, testTree.getRelPrec());
}

TEST_F(TreeArithmetic, SelfSubtraction) {
	{
		FunctionTree<D> subTree;
		subTree.add(1.0, *defTreeF, -1.0, *defTreeF);
		EXPECT_EQ(0.0, subTree.evalf(r1));
		EXPECT_EQ(0.0, subTree.evalf(r2));
		EXPECT_EQ(0.0, subTree.integrate());
		EXPECT_EQ(0.0, subTree.getSquareNorm());
	}
	{
		FunctionTree<D> subTree;
		subTree = *defTreeF - *defTreeF;
		EXPECT_EQ(0.0, subTree.evalf(r1));
		EXPECT_EQ(0.0, subTree.evalf(r2));
		EXPECT_EQ(0.0, subTree.integrate());
		EXPECT_EQ(0.0, subTree.getSquareNorm());
	}
	{
		FunctionTree<D> subTree;
		subTree = *defTreeF;
		subTree -= *defTreeF;
		EXPECT_EQ(0.0, subTree.evalf(r1));
		EXPECT_EQ(0.0, subTree.evalf(r2));
		EXPECT_EQ(0.0, subTree.integrate());
		EXPECT_EQ(0.0, subTree.getSquareNorm());
	}
}

TEST_F(TreeArithmetic, Addition) {
	GaussExp<D> refExp;
	refExp = *defFuncF + *defFuncG;

	FunctionTree<D> refTree;
	refTree.setRelPrec(refTree.getRelPrec() * 0.01);
	refTree.projectFunction(refExp);
	{
		FunctionTree<D> addTree;
		addTree.add(1.0, *defTreeF, 1.0, *defTreeG);
		testTrees(addTree, refTree);
	}
	{
		FunctionTree<D> addTree;
		addTree =  *defTreeF + *defTreeG;
		testTrees(addTree, refTree);
	}
	{
		FunctionTree<D> addTree;
		addTree = *defTreeF;
		addTree += *defTreeG;
		testTrees(addTree, refTree);
	}
}

TEST_F(TreeArithmetic, Multiplication) {
	GaussExp<D> refExp;
	refExp = *defFuncF * *defFuncG;

	FunctionTree<D> refTree;
	refTree.setRelPrec(refTree.getRelPrec() * 0.01);
	refTree.projectFunction(refExp);
	{
		FunctionTree<D> multTree;
		multTree.mult(1.0, *defTreeF, 1.0, *defTreeG);
		testTrees(multTree, refTree);
	}
	{
		FunctionTree<D> multTree;
		multTree =  *defTreeF * *defTreeG;
		testTrees(multTree, refTree);
	}
	{
		FunctionTree<D> multTree;
		multTree = *defTreeF;
		multTree *= *defTreeG;
		testTrees(multTree, refTree);
	}
}

TEST_F(TreeArithmetic, OperatorNesting) {
	GaussExp<D> refExp;
	refExp = *defFuncF * *defFuncG;

	FunctionTree<D> refTree;
	refTree.setRelPrec(refTree.getRelPrec() * 0.01);
	refTree.projectFunction(refExp);

	{
		FunctionTree<D> nestTree;
//		nestTree = 2.0*(*defTreeF)*(*defTreeG)-(*defTreeF)*(*defTreeG);

		FunctionTree<D> foo1;
		FunctionTree<D> foo2;

		foo1 = (*defTreeF) * (*defTreeG);
		foo1 *= 2.0;
		foo2 = (*defTreeF) * (*defTreeG);
		nestTree = foo1 - foo2;

		testTrees(nestTree, refTree);
	}

}

GaussExp<D> *TreeArithmetic::defFuncF = 0;
GaussExp<D> *TreeArithmetic::defFuncG = 0;
FunctionTree<D> *TreeArithmetic::defTreeF = 0;
FunctionTree<D> *TreeArithmetic::defTreeG = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "TreeArithmetic");
	FunctionTree<D>::setDefaultOrder(6);
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}

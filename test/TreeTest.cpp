#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "GaussConvolution.h"
#include "PoissonKernel.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "LebesgueIterator.h"

#define D 3

using namespace std;

class Environment: public ::testing::Environment {
public:
	virtual ~Environment() {}
	virtual void SetUp() {
	}
	virtual void TearDown() {
	}
};

class FuncTree: public ::testing::Test {
public:
	static void SetUpTestCase() {
		cout << scientific << setprecision(14);

		double r1[3] = {0.1, 0.2, 0.3};

		GaussExp<3> gtmp = MRUtils::genGaussExp<3>(3, -0.2, 0.1, 100.0);
		gfunc = new GaussFunc<3>(1234.0, 0.0, r1);
		gexp = new GaussExp<3>(gtmp);
		refTree = new FunctionTree<3>();
		refTree->projectFunction(*gfunc);
		refTree->setName("refTree");
	}
	static void TearDownTestCase() {
		delete gfunc;
		delete gexp;
		delete refTree;
	}

	virtual void SetUp() {
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

	static FunctionTree<3> *refTree;
	static GaussFunc<3> *gfunc;
	static GaussExp<3> *gexp;
};

TEST_F(FuncTree, EmptyTreeIsZero) {
	FunctionTree<3> tree;
	EXPECT_EQ(0.0, tree.evalf(r1));
	EXPECT_EQ(0.0, tree.evalf(r2));
}

TEST_F(FuncTree, ProjectGaussFunc) {
	double v1 = FuncTree::gfunc->evalf(r1);
	double v2 = FuncTree::refTree->evalf(r1);
	double dv = fabs(v1 - v2);

	EXPECT_LT(dv, FuncTree::refTree->getRelPrec());

}

TEST_F(FuncTree, ProjectGaussExp) {
	FuncTree::refTree->clear();
	FuncTree::refTree->projectFunction(*FuncTree::gexp);

	double v1 = FuncTree::gexp->evalf(r1);
	double v2 = FuncTree::refTree->evalf(r1);
	double dv = fabs(v1 - v2);

	EXPECT_LT(dv, FuncTree::refTree->getRelPrec()*0.1);
}

TEST_F(FuncTree, SaveTree) {
	FuncTree::refTree->saveTree("refTree");;
	println(0, *FuncTree::refTree);
}

TEST_F(FuncTree, LoadTree) {
	FunctionTree<3> fTree;
	fTree.loadTree("refTree");
	println(0, fTree);
	if (not fTree.isScattered()) {
		int diff = FuncTree::refTree->diffTree(fTree);
		EXPECT_EQ(0, diff);
	}
}

TEST_F(FuncTree, PurgeGenNodes) {
	FunctionTree<3> fTree;
	fTree.setAutoClean(true);
	fTree.projectFunction(*gexp);

	int totNodes[4] = {0,0,0,0};
	int projNodes[4] = {0,0,0,0};
	int endNodes[4] = {0,0,0,0};
	int genNodes[4] = {0,0,0,0};
	int allocGenNodes[4] = {0,0,0,0};

	{
		LebesgueIterator<3> it(&fTree);
		while (it.next()) { totNodes[0]++; }
	}
	projNodes[0] = fTree.getNNodes();
	endNodes[0] = fTree.getNEndNodes();
	genNodes[0] = fTree.getNGenNodes();
	allocGenNodes[0] = fTree.getNAllocGenNodes();

	EXPECT_LE(endNodes[0], projNodes[0]);
	EXPECT_EQ(0, totNodes[0] - projNodes[0]);
	EXPECT_EQ(0, genNodes[0]);
	EXPECT_EQ(0, allocGenNodes[0]);

	int n = fTree.getDepth();
	int l[3] = {0, 0, 0};
	NodeIndex<3> idx(n+3, l);

	MWNode<3> &node = fTree.getNode(idx, false);
	ASSERT_EQ(1, node.isGenNode());
	ASSERT_EQ(1, node.isAllocated());
	ASSERT_EQ(1, node.hasCoefs());

	{
		LebesgueIterator<3> it(&fTree);
		while (it.next()) { totNodes[1]++; }
	}
	projNodes[1] = fTree.getNNodes();
	endNodes[1] = fTree.getNEndNodes();
	genNodes[1] = fTree.getNGenNodes();
	allocGenNodes[1] = fTree.getNAllocGenNodes();

	EXPECT_EQ(0, endNodes[1] - endNodes[0]);
	EXPECT_EQ(0, projNodes[1] - projNodes[0]);
	EXPECT_LE(0, totNodes[1] - totNodes[0]);
	EXPECT_LE(0, genNodes[1]);
	EXPECT_LE(0, allocGenNodes[1]);
	EXPECT_LE(0, genNodes[1] - allocGenNodes[1]);

	fTree.clearGenNodes();

	{
		LebesgueIterator<3> it(&fTree);
		while (it.next()) { totNodes[2]++; }
	}
	projNodes[2] = fTree.getNNodes();
	endNodes[2] = fTree.getNEndNodes();
	genNodes[2] = fTree.getNGenNodes();
	allocGenNodes[2] = fTree.getNAllocGenNodes();

	EXPECT_EQ(0, endNodes[2] - endNodes[0]);
	EXPECT_EQ(0, projNodes[2] - projNodes[0]);
	EXPECT_EQ(0, totNodes[2] - totNodes[1]);
	EXPECT_LE(0, genNodes[2]);
	EXPECT_EQ(0, genNodes[2] - genNodes[1]);
	EXPECT_EQ(0, allocGenNodes[2]);

	fTree.purgeGenNodes();

	{
		LebesgueIterator<3> it(&fTree);
		while (it.next()) { totNodes[3]++; }
	}
	projNodes[3] = fTree.getNNodes();
	endNodes[3] = fTree.getNEndNodes();
	genNodes[3] = fTree.getNGenNodes();
	allocGenNodes[3] = fTree.getNAllocGenNodes();

	EXPECT_EQ(0, endNodes[3] - endNodes[0]);
	EXPECT_EQ(0, projNodes[3] - projNodes[0]);
	EXPECT_EQ(0, totNodes[3] - totNodes[0]);
	EXPECT_EQ(0, genNodes[3]);
	EXPECT_EQ(0, allocGenNodes[3]);
}

FunctionTree<3> *FuncTree::refTree = 0;
GaussFunc<3> *FuncTree::gfunc = 0;
GaussExp<3> *FuncTree::gexp = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
//	::testing::AddGlobalTestEnvironment(new Environment);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "FuncTreeTest");
	return RUN_ALL_TESTS();
}


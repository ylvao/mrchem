#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "GaussConvolution.h"
#include "PoissonKernel.h"
#include "OperatorTree.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"

using namespace std;

class Environment: public ::testing::Environment {
public:
	virtual ~Environment() {}
	virtual void SetUp() {
	}
	virtual void TearDown() {
	}
};

class Serialization: public ::testing::Test {
public:
	static void SetUpTestCase() {
		cout << scientific << setprecision(14);

		double r1[3] = {0.1, 0.2, 0.3};

		gfunc = new GaussFunc<3>(1234.0, 1.0, r1);
		refTree = new FunctionTree<3>();
		refTree->projectFunction(*gfunc);
		refTree->setName("refTree");

		PoissonKernel kernel(BoundingBox<3>::getWorldBox());
		refPoisson = new GaussConvolution<3>(kernel);
	}
	static void TearDownTestCase() {
		delete gfunc;
		delete refTree;
		delete refPoisson;
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
	static GaussConvolution<3> *refPoisson;
};

TEST_F(Serialization, SaveFunctionTree) {
	Serialization::refTree->saveTree("refTree");;
	println(0, *Serialization::refTree);
}

TEST_F(Serialization, LoadFunctionTree) {
	FunctionTree<3> fTree;
	fTree.loadTree("refTree");
	println(0, fTree);
	if (not fTree.isScattered()) {
		int diff = Serialization::refTree->diffTree(fTree);
		if (diff > 0) {
			println(1, "!!! Trees differ: " << diff << endl);
		} else {
			println(1, "*** Trees are identical. Hoorah!" << endl);
		}
		EXPECT_EQ(diff, 0);
	}
}

TEST_F(Serialization, SaveOperator) {
	refPoisson->saveOperator("refPoisson");
}

TEST_F(Serialization, LoadOperator) {
	GaussConvolution<3> loadedPoisson;
	loadedPoisson.loadOperator("refPoisson");

	FunctionTree<3> rgTree;
	rgTree = Serialization::refPoisson->apply(*Serialization::refTree);

	FunctionTree<3> lgTree;
	lgTree = loadedPoisson.apply(*Serialization::refTree);

	if (not rgTree.isScattered()) {
		int diff = rgTree.diffTree(lgTree);
		EXPECT_EQ(0, diff);
	}

	println(0, "------------------ ref. gTree --------------");
	println(0, rgTree);
	println(0, "------------------ load gTree --------------");
	println(0, lgTree);
	println(0, "---------------------------------------------");
	double rint = rgTree.integrate();
	double lint = lgTree.integrate();
	println(0, "ref.   gTree integral: " << rint);
	println(0, "loaded gTree integral: " << lint);
	EXPECT_LT(fabs(rint-lint), 10e-9);
}

FunctionTree<3> *Serialization::refTree = 0;
GaussConvolution<3> *Serialization::refPoisson = 0;
GaussFunc<3> *Serialization::gfunc = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
//	::testing::AddGlobalTestEnvironment(new Environment);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "");
	return RUN_ALL_TESTS();
}

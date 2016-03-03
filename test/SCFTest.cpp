#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "CoulombOperator.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "HartreeFock.h"
#include "XCFun.h"
#include "DFT.h"
#include "Orbital.h"
#include "Molecule.h"

CoulombOperator *defaultCoulomb;
using namespace std;

class SCFTest: public ::testing::Test {
public:
	static void SetUpTestCase() {
		defaultCoulomb = new CoulombOperator();
	}
	static void TearDownTestCase() {
		delete defaultCoulomb;
	}

	virtual void SetUp() {
		SET_PRINT_PRECISION(15);
		cout << scientific << setprecision(15);
		rank = 0;
		mpi::communicator world;
		rank = world.rank();

		string basis_file = "he.bas";
		string dens_file = Input.get<string>("Molecule.density");
		string overlap_file = Input.get<string>("Molecule.overlap");

		helium = new Molecule(basis_file, dens_file, overlap_file);
		helium->initializeNuclearPotential(Molecule::HarrisonTypePotential);
	}

	virtual void TearDown() {
		delete helium;
		helium = 0;
	}

	int rank;

	static Molecule *helium;

};

TEST_F(SCFTest, OneElectron) {
	helium->initializeOrbitals(2, 1);

	HartreeFock hf;
	hf.setPrec(1.0e-1);
	hf.optimizeOrbital(*helium);

	double energy;

	energy = helium->getOrbitals().getOrbital(0).getEnergy();
	println(0, "Energy: " << energy);
	EXPECT_LT(fabs(energy + 2.0), 0.01);

	energy = helium->getOrbitals().getOrbital(1).getEnergy();
	println(0, "Energy: " << energy);
	EXPECT_LT(fabs(energy + 0.5), 0.01);
}

TEST_F(SCFTest, TwoElectronDFT) {
	helium->initializeOrbitals(1, 2);

	XCFun xcfun;
	DFT dft(xcfun);
	dft.setPrec(1.0e-2);
	dft.optimizeOrbital(*helium);

	double energy = helium->getOrbitals().getOrbital(0).getEnergy();
	println(0, "Energy: " << energy);
	EXPECT_LT(fabs(energy + 0.57), 0.01);
}

TEST_F(SCFTest, TwoElectronHF) {
	helium->initializeOrbitals(1, 2);

	HartreeFock hf;
	hf.setPrec(1.0e-2);
	hf.optimizeOrbital(*helium);

	double energy = helium->getOrbitals().getOrbital(0).getEnergy();
	println(0, "Energy: " << energy);
	EXPECT_LT(fabs(energy + 0.91), 0.01);
}

Molecule *SCFTest::helium = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "SCFTest");
	SET_PRINT_PRECISION(15);
	return RUN_ALL_TESTS();
}

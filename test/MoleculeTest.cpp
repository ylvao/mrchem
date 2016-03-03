#include "GaussFunc.h"
#include "GaussExp.h"
#include "Plot.h"
#include "GaussConvolution.h"
#include "PoissonKernel.h"
#include "Molecule.h"
#include "mwtest.h"
#include "MRUtils.h"
#include "NuclearGaussPotential.h"
#include "NuclearPotential.h"

using namespace std;

CoulombOperator *defaultCoulomb;

class MoleculeTest: public ::testing::Test {
public:
	static void SetUpTestCase() {
		cout << scientific << setprecision(14);

		string basis_file = Input.get<string>("Molecule.basis");
		string dens_file = Input.get<string>("Molecule.density");
		string overlap_file = Input.get<string>("Molecule.overlap");

		harrisonMolecule = new Molecule(basis_file, dens_file, overlap_file);
		gaussMolecule = new Molecule(basis_file, dens_file, overlap_file);

		defaultCoulomb = new CoulombOperator();
	}
	static void TearDownTestCase() {
		delete harrisonMolecule;
		delete gaussMolecule;
		delete defaultCoulomb;
	}

	virtual void SetUp() {
		rank = 0;
		mpi::communicator world;
		rank = world.rank();
	}

	virtual void TearDown() {
	}

	int rank;

	static Molecule *harrisonMolecule;
	static Molecule *gaussMolecule;
};

TEST_F(MoleculeTest, GaussPotential) {
	double nucExp = 10;

	gaussMolecule->setNuclearExponent(nucExp);
	gaussMolecule->initializeNuclearPotential(Molecule::GaussTypePotential);

	NuclearGaussPotential &nucGaussPot = static_cast<NuclearGaussPotential &>
				(gaussMolecule->getNuclearPotential());

	Density nucCharge;
	nucCharge.projectFunction(nucGaussPot.getNucExpansion());
	nucCharge.calcNumberDensity();
	double numDens = nucCharge.getNumberDensity();
	double prec = nucCharge.getRelPrec();

	println(0, "Nuclear charge: " << numDens);

	EXPECT_LT(fabs(numDens - 10), prec);
}

TEST_F(MoleculeTest, HarrisonPotential) {
	harrisonMolecule->initializeNuclearPotential(Molecule::HarrisonTypePotential);

	NuclearPotential &harrisonPot = harrisonMolecule->getNuclearPotential();
	NuclearPotential &gaussPot = gaussMolecule->getNuclearPotential();

	double harrisonInt = harrisonPot.integrate();
	double gaussInt = gaussPot.integrate();

	double relDiff = fabs(harrisonInt - gaussInt)/gaussPot.getSquareNorm();

	EXPECT_LT(relDiff, 1.0);
}

TEST_F(MoleculeTest, ElectronDensity) {
	harrisonMolecule->calcAODensity();

	Density elDens = harrisonMolecule->getElectronDensity();
	elDens.calcNumberDensity();

	double prec = elDens.getRelPrec();
	double numDens = elDens.getNumberDensity();
	println(0, "Number density: " << numDens);

	EXPECT_LT(fabs(numDens - 10), prec);
}

TEST_F(MoleculeTest, ElectronicPotential) {
	Density &elDens = harrisonMolecule->getElectronDensity();

	CoulombPotential elPot;
	elPot.calcPotential(elDens);

	const double refEnergy = 93.72165;
	double coulombEnergy = elPot.calcEnergy(elDens);
	println(0, "Coulomb energy: " << coulombEnergy);

	EXPECT_LT(fabs(coulombEnergy - refEnergy), 0.1);
}

Molecule *MoleculeTest::harrisonMolecule = 0;
Molecule *MoleculeTest::gaussMolecule = 0;

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	mpi::environment env(argc, argv);
	MREnv::initializeMRCPP(argc, argv, "MoleculeTest");
	return RUN_ALL_TESTS();
}




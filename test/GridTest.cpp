#include "mwtest.h"
#include "BoundingBox.h"
#include "Molecule.h"
#include "AtomicElement.h"
#include "PeriodicTable.h"
#include "MolecularGridGenerator.h"
#include "MRGrid.h"

using namespace std;
using namespace Eigen;

class GridTest: public ::testing::Test {
public:
    static void SetUpTestCase() {
	SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
    }

    virtual void TearDown() {
    }
};

TEST_F(GridTest, Molecule) {
    PeriodicTable pt;
    const AtomicElement &O = pt.getAtomicElement("O");
    const AtomicElement &H = pt.getAtomicElement("H");

    double r[3];
    Molecule water;
    Atom oxygen(O);
    Atom hydrogen(H);

    r[0] = 0.0000;
    r[1] = 0.0000;
    r[2] = 0.0000;
    oxygen.setCoord(r);
    water.addAtom(oxygen);

    r[0] = 1.4375;
    r[1] = 0.0000;
    r[2] = 1.1500;
    hydrogen.setCoord(r);
    water.addAtom(hydrogen);

    r[0] =-1.4375;
    r[1] = 0.0000;
    r[2] = 1.1500;
    hydrogen.setCoord(r);
    water.addAtom(hydrogen);

    const double *h_coord = water.getAtom(2).getCoord();

    EXPECT_EQ(water.getNAtoms(), 3);
    EXPECT_LT(fabs(h_coord[2]-1.15), MachineZero);

    println(0, endl);
    println(0, water);
    println(0, endl);
}

TEST_F(GridTest, UniformGrid_1D) {
    int order = 5;
    int uniform = 3;

    int rootScale = -2;
    int nBoxes[1] = {3};
    double origin[1] = {6.0};
    BoundingBox<1> box_1D(rootScale, nBoxes, origin);
    println(0, box_1D);

    MRGrid<1> grid_1D(order, &box_1D);
    const double *lb = grid_1D.getLowerBounds();
    const double *ub = grid_1D.getUpperBounds();
    EXPECT_EQ(grid_1D.getNNodes(), 3);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 

/*
    GridGenerator<1> gen_1D(uniform);
    gen_1D.generateGrid(grid_1D);

    EXPECT_EQ(grid_1D.getNNodes(), 27);
    EXPECT_EQ(grid_1D.getNQuadraturePoints(), 324);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
*/
}

TEST_F(GridTest, UniformGrid_2D) {
    int order = 3;
    int uniform = 2;

    int rootScale = -2;
    int nBoxes[2] = {3,2};
    double origin[2] = {6.0,4.0};
    BoundingBox<2> box_2D(rootScale, nBoxes, origin);
    println(0, box_2D);
/*
    MRGrid<2> grid_2D(order, box_2D);
    const double *lb = grid_2D.getLowerBounds();
    const double *ub = grid_2D.getUpperBounds();
    const double *l = grid_2D.getLength();
    EXPECT_EQ(grid_2D.getNNodes(), 6);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 

    GridGenerator<2> gen_2D(uniform);
    gen_2D.generateGrid(grid_2D);
    EXPECT_EQ(grid_2D.getNNodes(), 96);
    EXPECT_EQ(grid_2D.getNQuadraturePoints(), 3456);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
*/
}

TEST_F(GridTest, UniformGrid_3D) {
    int order = 5;
    int uniform = 1;

    int rootScale = -2;
    int nBoxes[3] = {3,2,1};
    double origin[3] = {6.0,4.0,2.0};
    BoundingBox<3> box_3D(rootScale, nBoxes, origin);
    println(0, box_3D);
/*
    MRGrid<3> grid_3D(order, box_3D);
    const double *lb = grid_2D.getLowerBounds();
    const double *ub = grid_2D.getUpperBounds();
    const double *l = grid_2D.getLength();
    EXPECT_EQ(grid_2D.getNNodes(), 6);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(lb[2] + 2.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[2] - 2.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
    EXPECT_LT(fabs(l[2] -  4.0), MachineZero); 

    GridGenerator<3> gen_3D(uniform);
    gen_3D.generateGrid(grid_3D);
    EXPECT_EQ(grid_2D.getNNodes(), 48);
    EXPECT_EQ(grid_2D.getNQuadraturePoints(), 1728);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(lb[2] + 2.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[2] - 2.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
    EXPECT_LT(fabs(l[2] -  4.0), MachineZero); 
*/
}

TEST_F(GridTest, MolecularGrid) {
    PeriodicTable pt;
    const AtomicElement &H = pt.getAtomicElement("H");

    double r1[3] = {0.0, 0.0, -0.7};
    double r2[3] = {0.0, 0.0,  0.7};
    
    Atom h1(H, r1);
    Atom h2(H, r2);

    Molecule H2;
    H2.addAtom(h1);
    H2.addAtom(h2);

    EXPECT_EQ(H2.getNAtoms(), 2);

    println(0, endl);
    println(0, H2);
    println(0, endl);

    MolecularGridGenerator generator;
    generator.setUniform(3);
    generator.setDepth(5);
    generator.setWidth(5);

    //MRGrid<3> grid;
    //generator.generateGrid(grid, H2);

    //VectorXd weights = grid.getQuadratureWeights();
    //MatrixXd points = grid.getQuadraturePoints();

    //println(0, grid);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::AddGlobalTestEnvironment(new Environment);
    mpi::environment env(argc, argv);
    //MREnv::initializeMRCPP(argc, argv, "FuncTreeTest");
    return RUN_ALL_TESTS();
}



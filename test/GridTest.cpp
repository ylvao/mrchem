#include "mwtest.h"
#include "NodeBox.h"
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

    //println(0, endl);
    //println(0, water);
    //println(0, endl);
}

TEST_F(GridTest, UniformGrid_1D) {
    int rootScale = -2;
    int nBoxes[1] = {3};
    double origin[1] = {6.0};
    NodeBox<1> box_1D(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<1> grid_1D(order, &box_1D);
    const double *lb = grid_1D.getLowerBounds();
    const double *ub = grid_1D.getUpperBounds();
    const double *l = grid_1D.getRootBox().getBoxLength();
    EXPECT_EQ(grid_1D.getNNodes(), 3);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 

    int uniform = 1;
    GridGenerator<1> gen_1D;
    gen_1D.setUniformScale(uniform);
    gen_1D.generateGrid(grid_1D);
    EXPECT_EQ(grid_1D.getNNodes(), 45);
    EXPECT_EQ(grid_1D.countLeafNodes(), 24);
    EXPECT_EQ(grid_1D.countQuadPoints(), 288);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 

    uniform = 5;
    gen_1D.setUniformScale(uniform);
    gen_1D.generateGrid(grid_1D);
    //println(0, grid_1D);

    VectorXd weights;
    MatrixXd points;

    grid_1D.getQuadWeights(weights);
    grid_1D.getQuadPoints(points);

    double alpha = 1.0;
    double coef = sqrt(alpha/pi);
    double pos = 0.0;

    int nPoints = points.rows();
    VectorXd values = VectorXd::Zero(nPoints);
    for (int i = 0; i < nPoints; i++) {
	double r = pos - points(i,0);
	values(i) = coef*exp(-alpha*r*r);
    }

    double integral = weights.dot(values);
    //println(0, "Integral " << integral);

}

TEST_F(GridTest, UniformGrid_2D) {
    int rootScale = -2;
    int nBoxes[2] = {3,2};
    double origin[2] = {6.0,4.0};
    NodeBox<2> box_2D(rootScale, nBoxes, origin);

    int order = 3;
    MRGrid<2> grid_2D(order, &box_2D);
    const double *lb = grid_2D.getLowerBounds();
    const double *ub = grid_2D.getUpperBounds();
    const double *l = grid_2D.getRootBox().getBoxLength();
    EXPECT_EQ(grid_2D.getNNodes(), 6);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 

    int uniform = 0;
    GridGenerator<2> gen_2D;
    gen_2D.setUniformScale(uniform);
    gen_2D.generateGrid(grid_2D);
    EXPECT_EQ(grid_2D.getNNodes(), 126);
    EXPECT_EQ(grid_2D.countLeafNodes(), 96);
    EXPECT_EQ(grid_2D.countQuadPoints(), 6144);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
}

TEST_F(GridTest, UniformGrid_3D) {
    int rootScale = -2;
    int nBoxes[3] = {3,2,1};
    double origin[3] = {6.0,4.0,2.0};
    NodeBox<3> box_3D(rootScale, nBoxes, origin);

    int order = 5;
    MRGrid<3> grid_3D(order, &box_3D);
    const double *lb = grid_3D.getLowerBounds();
    const double *ub = grid_3D.getUpperBounds();
    const double *l = grid_3D.getRootBox().getBoxLength();
    EXPECT_EQ(grid_3D.getNNodes(), 6);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(lb[2] + 2.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[2] - 2.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
    EXPECT_LT(fabs(l[2] -  4.0), MachineZero); 

    int uniform = -1;
    GridGenerator<3> gen_3D;
    gen_3D.setUniformScale(uniform);
    gen_3D.generateGrid(grid_3D);
    EXPECT_EQ(grid_3D.getNNodes(), 54);
    EXPECT_EQ(grid_3D.countLeafNodes(), 48);
    EXPECT_EQ(grid_3D.countQuadPoints(), 82944);
    EXPECT_LT(fabs(lb[0] + 6.0), MachineZero); 
    EXPECT_LT(fabs(lb[1] + 4.0), MachineZero); 
    EXPECT_LT(fabs(lb[2] + 2.0), MachineZero); 
    EXPECT_LT(fabs(ub[0] - 6.0), MachineZero); 
    EXPECT_LT(fabs(ub[1] - 4.0), MachineZero); 
    EXPECT_LT(fabs(ub[2] - 2.0), MachineZero); 
    EXPECT_LT(fabs(l[0] - 12.0), MachineZero); 
    EXPECT_LT(fabs(l[1] -  8.0), MachineZero); 
    EXPECT_LT(fabs(l[2] -  4.0), MachineZero); 
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

    //println(0, endl);
    //println(0, H2);
    //println(0, endl);

    int rootScale = -3;
    int nBoxes[3] = {1,1,1};
    double origin[3] = {4.0,4.0,4.0};
    NodeBox<3> world(rootScale, nBoxes, origin);
    //println(0, world);

    MolecularGridGenerator generator;
    generator.setUniformScale(-1);
    generator.setAmplitude(5);
    generator.setWidth(5);
    generator.setNuclearDependence(0);

    int order = 3;
    MRGrid<3> grid(order, &world);
    generator.generateGrid(grid, H2);
    //println(0, grid);
}

TEST_F(GridTest, GridIntegral) {
    PeriodicTable pt;
    const AtomicElement &H = pt.getAtomicElement("H");
    double pos[3] = {pi, pi, pi};
    
    Atom h(H, pos);
    //println(0, endl);
    //println(0, h);
    //println(0, endl);

    int rootScale = -3;
    int nBoxes[3] = {1,1,1};
    double origin[3] = {0.0,0.0,0.0};
    NodeBox<3> world(rootScale, nBoxes, origin);
    //println(0, endl);
    //println(0, world);
    //println(0, endl);

    MolecularGridGenerator generator;
    generator.setUniformScale(-1);
    generator.setAmplitude(3);
    generator.setWidth(1);
    generator.setNuclearDependence(0);

    int order = 3;
    MRGrid<3> grid(order, &world);
    generator.generateGrid(grid, h);
    //println(0, grid);

    VectorXd weights;
    MatrixXd points;

    grid.getQuadWeights(weights);
    grid.getQuadPoints(points);

    double alpha = 10.0;
    double coef = pow(alpha/pi, 3.0/2.0);

    int nPoints = points.rows();
    VectorXd values = VectorXd::Zero(nPoints);
    for (int i = 0; i < nPoints; i++) {
	values(i) = coef;
	for (int d = 0; d < 3; d++) {
	    double r = pos[d] - points(i,d);
	    values(i) *= exp(-alpha*r*r);
	}
    }

    double integral = weights.dot(values);
    EXPECT_LT(integral - 1.0, 1.0e-5);
    //println(0, endl);
    //println(0, "Integral " << integral);
    //println(0, endl);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    //::testing::AddGlobalTestEnvironment(new Environment);
    mpi::environment env(argc, argv);
    //MREnv::initializeMRCPP(argc, argv, "FuncTreeTest");
    return RUN_ALL_TESTS();
}



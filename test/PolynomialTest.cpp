#include "mwtest.h"
#include "Polynomial.h"

using namespace std;
using namespace Eigen;

class PolynomialTest: public ::testing::Test {
public:
    static void SetUpTestCase() {
	SET_PRINT_PRECISION(15);
    }
    static void TearDownTestCase() {
    }

    virtual void SetUp() {
	X = 0.3;
	A = -1.0;
	B = 1.0;

	F = 0;
	G = 0;

	initializeReferenceF();
	initializeReferenceG();
    }

    virtual void TearDown() {
	ASSERT_TRUE(F != 0);
	testBounds(F, &A, &B);
	testReferenceF(F);
	deletePolynomial(&F);
	ASSERT_TRUE(F == 0);

	ASSERT_TRUE(G != 0);
	testBounds(G, &A, &B);
	testReferenceG(G);
	deletePolynomial(&G);
	ASSERT_TRUE(G == 0);
    }

    void initializeReferenceF();
    void initializeReferenceG();

    void testReferenceF(const Polynomial *f);
    void testReferenceG(const Polynomial *g);

    void testBounds(const Polynomial *poly, const double *a, const double *b);
    void testPolynomial(const Polynomial *poly, const Polynomial *ref);
    void deletePolynomial(Polynomial **poly);

    double X;
    double A;
    double B;

    Polynomial *F;
    Polynomial *G;
};

void PolynomialTest::initializeReferenceF() {
    ASSERT_TRUE(F == 0);
    VectorXd coefs = VectorXd::Zero(4);
    coefs[0] = 0.0;
    coefs[1] = 1.0;
    coefs[2] = 1.0;
    coefs[3] = 0.0;
    F = new Polynomial(coefs, &A, &B);
}

void PolynomialTest::initializeReferenceG() {
    ASSERT_TRUE(G == 0);
    VectorXd coefs = VectorXd::Zero(2);
    coefs[0] = 0.0;
    coefs[1] = 1.0;
    G = new Polynomial(coefs, &A, &B);
}

void PolynomialTest::deletePolynomial(Polynomial **poly) {
    ASSERT_TRUE(poly != 0);
    ASSERT_TRUE(*poly != 0);
    delete *poly;
    *poly = 0;
}

void PolynomialTest::testBounds(const Polynomial *poly, const double *a, const double *b) {
    if (a == 0 or b == 0) {
	ASSERT_TRUE(a == b);
        EXPECT_FALSE(poly->isBounded());
	EXPECT_TRUE(poly->getLowerBounds() == 0);
	EXPECT_TRUE(poly->getUpperBounds() == 0);
    } else {
	ASSERT_TRUE(a != 0);
	ASSERT_TRUE(b != 0);
        EXPECT_TRUE(poly->isBounded());
	EXPECT_TRUE(poly->getLowerBounds() != 0);
	EXPECT_TRUE(poly->getUpperBounds() != 0);
	EXPECT_LT(fabs(poly->getLowerBounds()[0] - *a), MachineZero);
	EXPECT_LT(fabs(poly->getUpperBounds()[0] - *b), MachineZero);
    }
}

void PolynomialTest::testReferenceF(const Polynomial *f) {
    ASSERT_TRUE(f != 0);

    EXPECT_EQ(4, f->size());
    EXPECT_EQ(2, f->getOrder());

    EXPECT_LT(f->getSquareNorm(), 0.0);
    EXPECT_LT(fabs(f->getDilation() - 1.0), MachineZero);
    EXPECT_LT(fabs(f->getTranslation()), MachineZero);

    const double val = X + X*X;
    const double err = val - f->evalf(X);
    EXPECT_LT(fabs(err), MachineZero);
}

void PolynomialTest::testReferenceG(const Polynomial *g) {
    ASSERT_TRUE(g != 0);
}

void PolynomialTest::testPolynomial(const Polynomial *poly, const Polynomial *ref) {
    ASSERT_TRUE(poly != 0);
    ASSERT_TRUE(ref != 0);

    EXPECT_EQ(ref->size(), poly->size());
    EXPECT_EQ(ref->getOrder(), poly->getOrder());

    const double normDiff = poly->getSquareNorm() - ref->getSquareNorm();
    const double dilDiff = poly->getDilation() - ref->getDilation();
    const double translDiff = poly->getTranslation() - ref->getTranslation();
    const double valDiff = poly->evalf(X) - ref->evalf(X);

    EXPECT_LT(fabs(normDiff), MachineZero);
    EXPECT_LT(fabs(dilDiff), MachineZero);
    EXPECT_LT(fabs(translDiff), MachineZero);
    EXPECT_LT(fabs(valDiff), MachineZero);
}

TEST_F(PolynomialTest, Initialization) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
}

TEST_F(PolynomialTest, DefaultConstructor) {
    Polynomial *f = new Polynomial;
    testBounds(f, 0, 0);

    EXPECT_EQ(1, f->size());
    EXPECT_EQ(0, f->getOrder());

    EXPECT_LT(f->getSquareNorm(), 0.0);
    EXPECT_LT(fabs(f->getDilation() - 1.0), MachineZero);
    EXPECT_LT(fabs(f->getTranslation()), MachineZero);
    EXPECT_LT(fabs(f->evalf(X)), MachineZero);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, CopyConstructor) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);

    testBounds(f, &A, &B);
    testReferenceF(f);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, AssignmentOperator) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial;
    *f = *F;

    testBounds(f, 0, 0);
    testReferenceF(f);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, SetCoefs) {
    VectorXd coefs = VectorXd::Zero(4);
    coefs[0] = 0.0;
    coefs[1] = 1.0;
    coefs[2] = 1.0;
    coefs[3] = 0.0;

    Polynomial *f = new Polynomial;
    f->setCoefs(coefs);

    testBounds(f, 0, 0);
    testReferenceF(f);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, OutOfBoundsEvalf) {
    ASSERT_TRUE(F != 0);
    EXPECT_LT(fabs(F->evalf(-2.000)), MachineZero);
    EXPECT_LT(fabs(F->evalf(-1.001)), MachineZero);
    EXPECT_LT(fabs(F->evalf(1.001)), MachineZero);
    EXPECT_LT(fabs(F->evalf(10.00)), MachineZero);
}

TEST_F(PolynomialTest, Rescale) {
    ASSERT_TRUE(F != 0);
    const double n = 2.0;
    const double l = 1.0;

    Polynomial *f = new Polynomial(*F);
    f->rescale(n, l);

    const double dilDiff = f->getDilation() - 2.0;
    const double translDiff = f->getTranslation() - 1.0;
    EXPECT_LT(fabs(dilDiff), MachineZero);
    EXPECT_LT(fabs(translDiff), MachineZero);

    const double aDiff = f->getLowerBound() - 0.00;
    const double bDiff = f->getUpperBound() - 1.00;
    EXPECT_LT(fabs(aDiff), MachineZero);
    EXPECT_LT(fabs(bDiff), MachineZero);

    const double val = (2.0*X - 1.0) + (2.0*X - 1.0)*(2.0*X - 1.0);
    const double valDiff = val - f->evalf(X);
    EXPECT_LT(fabs(valDiff), MachineZero);

    Polynomial *g = new Polynomial(*F);
    g->setDilation(n);
    g->setTranslation(l);

    testBounds(f, &A, &B);
    testBounds(g, &A, &B);
    testPolynomial(f, g);
    
    deletePolynomial(&f); 
    deletePolynomial(&g); 
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, ClearCoefs) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->clearCoefs();

    testBounds(f, &A, &B);
    EXPECT_EQ(1, f->size());
    EXPECT_EQ(0, f->getOrder());
    EXPECT_LT(fabs(f->evalf(X)), MachineZero);
    EXPECT_LT(f->getSquareNorm(), 0.0);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, SetZero) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->setZero();

    testBounds(f, &A, &B);
    EXPECT_EQ(4, f->size());
    EXPECT_EQ(0, f->getOrder());
    EXPECT_LT(fabs(f->evalf(X)), MachineZero);
    EXPECT_LT(f->getSquareNorm(), 0.0);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, MultConstInPlace) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    *f *= 3.0;

    Polynomial *g = new Polynomial(3, &A, &B);
    VectorXd &coefs = g->getCoefs();
    coefs[1] = 3.0;
    coefs[2] = 3.0;
    
    testBounds(f, &A, &B);
    testBounds(g, &A, &B);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, MultConst) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial;
    *f = (*F)*2.0;

    Polynomial *g = new Polynomial(3);
    VectorXd &coefs = g->getCoefs();
    coefs[1] = 2.0;
    coefs[2] = 2.0;

    testBounds(f, 0, 0);
    testBounds(g, 0, 0);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, MultInPlace) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
    Polynomial *g = new Polynomial(*G);
    *g *= *F;

    Polynomial *h = new Polynomial(3, &A, &B);
    VectorXd &coefs = h->getCoefs();
    coefs[2] = 1;
    coefs[3] = 1;

    testBounds(g, &A, &B);
    testBounds(h, &A, &B);
    testPolynomial(g, h);

    deletePolynomial(&g);
    deletePolynomial(&h);
    ASSERT_TRUE(g == 0);
    ASSERT_TRUE(h == 0);
}

TEST_F(PolynomialTest, Multiplication) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
    Polynomial *h = new Polynomial;
    *h = (*F)*(*G);

    Polynomial *fg = new Polynomial(3, &A, &B);
    VectorXd &coefs = fg->getCoefs();
    coefs[2] = 1.0;
    coefs[3] = 1.0;
    
    testBounds(h, 0, 0);
    testBounds(fg, &A, &B);
    testPolynomial(h, fg);

    deletePolynomial(&h);
    deletePolynomial(&fg);
    ASSERT_TRUE(h == 0);
    ASSERT_TRUE(fg == 0);
}

TEST_F(PolynomialTest, AddInPlace) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
    Polynomial *g = new Polynomial(*G);
    *g += *F;

    Polynomial *h = new Polynomial(2, &A, &B);
    VectorXd &coefs = h->getCoefs();
    coefs[1] = 2;
    coefs[2] = 1;

    testBounds(g, &A, &B);
    testBounds(h, &A, &B);
    testPolynomial(g, h);

    deletePolynomial(&g);
    deletePolynomial(&h);
    ASSERT_TRUE(g == 0);
    ASSERT_TRUE(h == 0);
}

TEST_F(PolynomialTest, Addition) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
    Polynomial *h = new Polynomial;
    *h = (*F)+(*G);

    Polynomial *fg = new Polynomial(2, &A, &B);
    VectorXd &coefs = fg->getCoefs();
    coefs[1] = 2.0;
    coefs[2] = 1.0;
    
    testBounds(h, 0, 0);
    testBounds(fg, &A, &B);
    testPolynomial(h, fg);

    deletePolynomial(&h);
    deletePolynomial(&fg);
    ASSERT_TRUE(h == 0);
    ASSERT_TRUE(fg == 0);
}

TEST_F(PolynomialTest, DerivativeInPlace) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->calcDerivativeInPlace();

    Polynomial *g = new Polynomial(1, &A, &B);
    VectorXd &coefs = g->getCoefs();
    coefs[0] = 1.0;
    coefs[1] = 2.0;

    testBounds(f, &A, &B);
    testBounds(g, &A, &B);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, Derivative) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial;
    *f = F->calcDerivative();

    Polynomial *g = new Polynomial(1);
    VectorXd &coefs = g->getCoefs();
    coefs[0] = 1.0;
    coefs[1] = 2.0;

    testBounds(f, 0, 0);
    testBounds(g, 0, 0);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, AntiDerivativeInPlace) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->calcAntiDerivativeInPlace();

    Polynomial *g = new Polynomial(3, &A, &B);
    VectorXd &coefs = g->getCoefs();
    coefs[2] = 1.0/2.0;
    coefs[3] = 1.0/3.0;

    testBounds(f, &A, &B);
    testBounds(g, &A, &B);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, AntiDerivative) {
    ASSERT_TRUE(F != 0);
    Polynomial *f = new Polynomial;
    *f = F->calcAntiDerivative();

    Polynomial *g = new Polynomial(3);
    VectorXd &coefs = g->getCoefs();
    coefs[2] = 1.0/2.0;
    coefs[3] = 1.0/3.0;

    testBounds(f, 0, 0);
    testBounds(g, 0, 0);
    testPolynomial(f, g);

    deletePolynomial(&f);
    deletePolynomial(&g);
    ASSERT_TRUE(f == 0);
    ASSERT_TRUE(g == 0);
}

TEST_F(PolynomialTest, Integrate) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);

    const double a = 0.0;
    const double b = 1.0;

    const double fullErr = F->integrate() - 2.0/3.0;
    const double partErr = G->integrate(&a, &b) - 1.0/2.0;
    
    EXPECT_LT(fabs(fullErr), MachineZero);
    EXPECT_LT(fabs(partErr), MachineZero);
}

TEST_F(PolynomialTest, InnerProduct) {
    ASSERT_TRUE(F != 0);
    ASSERT_TRUE(G != 0);
    const double innerProd = F->innerProduct(*G);
    const double err = innerProd - 2.0/3.0;
    EXPECT_LT(fabs(err), MachineZero);
}

TEST_F(PolynomialTest, CalcSquareNorm) {
    EXPECT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->calcSquareNorm();

    const double normDiff = f->getSquareNorm() - 16.0/15.0;
    EXPECT_LT(fabs(normDiff), MachineZero);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

TEST_F(PolynomialTest, Normalize) {
    EXPECT_TRUE(F != 0);
    Polynomial *f = new Polynomial(*F);
    f->normalize();

    const double normDiff = f->getSquareNorm() - 1.0;
    EXPECT_LT(fabs(normDiff), MachineZero);

    deletePolynomial(&f);
    ASSERT_TRUE(f == 0);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
//  ::testing::AddGlobalTestEnvironment(new Environment);
//  MREnv::initializeMRCPP(argc, argv, "FuncTreeTest");
    mpi::environment env(argc, argv);
    return RUN_ALL_TESTS();
}



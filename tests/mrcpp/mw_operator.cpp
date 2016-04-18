#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "IdentityOperator.h"

namespace mw_projector {

template<int D> void testIdentityOperator();

SCENARIO("Applying identity operator", "[mw_operator], [tree_builder], [trees]") { 
    GIVEN("a MW representation of a Gaussian of unit charge in 1D") {
       testIdentityOperator<1>();
    }
    GIVEN("a MW representation of a Gaussian of unit charge in 2D") {
        testIdentityOperator<2>();
    }
    GIVEN("a MW representation of a Gaussian of unit charge in 3D") {
        testIdentityOperator<3>();
    }
}

template<int D> void testIdentityOperator() {
    GaussFunc<D> *func = 0;
    initialize(&func);
    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);
	MWProjector<D> Q(*mra,1.0e-5);
	FunctionTree<D> *f = Q(*func);

    WHEN("the identity is applied on the default grid") {
		IdentityOperator<D> I(*mra,1.0e-3);
		GridGenerator<D> G(*mra);
		FunctionTree<D> *g = G();
		I(*g,*f);
		println(0,f->getNNodes());
		println(0,g->getNNodes());
        THEN("it integrates to approximately one") {
            REQUIRE( f->integrate() == Approx(g->integrate()) );			
		}
	}
}
	/*
            REQUIRE( tree->integrate() == Approx(1.0).epsilon(1.0e+1) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree->getSquareNorm();
            REQUIRE( tree->dot(*tree) == Approx(norm) );
        }
        delete tree;
    }
    WHEN("the function is projected on an adapted grid") {
        GridGenerator<D> G(*mra);
        MWProjector<D> Q(*mra);
        FunctionTree<D> *tree = G(*func);
        Q(*tree, *func);
        THEN("it integrates to approximately one") {
            REQUIRE( tree->integrate() == Approx(1.0).epsilon(1.0e-3) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = tree->getSquareNorm();
            REQUIRE( tree->dot(*tree) == Approx(norm) );
        }
        delete tree;
    }
    WHEN("the function is projected with guaranteed precision") {
        const double prec = 1.0e-4;
        MWProjector<D> Q_adap(*mra, prec);
        FunctionTree<D> *f_tree = Q_adap(*func);
        THEN("it integrates to approximately one") {
            REQUIRE( f_tree->integrate() == Approx(1.0).epsilon(1.0e-8) );
        }
        THEN("the dot product with itself is equal to its squared norm") {
            const double norm = f_tree->getSquareNorm();
            REQUIRE( f_tree->dot(*f_tree) == Approx(norm) );
        }
        AND_WHEN("the function is projected on an identical grid") {
            CopyAdaptor<D> c_adap(*f_tree);
            MWProjector<D> Q_copy(*mra, c_adap);
            FunctionTree<D> *g_tree = Q_copy(*func);
            THEN("it integrates to the same value") {
                const double charge = f_tree->integrate();
                REQUIRE( g_tree->integrate() == Approx(charge) );
            }
            THEN("the dot product with the original is equal to their squared norm") {
                const double norm = f_tree->getSquareNorm();
                REQUIRE( g_tree->dot(*f_tree) == Approx(norm) );
            }
            delete g_tree;
        }
        delete f_tree;
    }
    finalize(&mra);
    finalize(&func);
}
	*/
} // namespace

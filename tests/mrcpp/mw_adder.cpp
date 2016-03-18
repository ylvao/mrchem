#include "catch.hpp"

#include "factory_functions.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "MWAdder.h"
#include "WaveletAdaptor.h"

namespace mw_adder {

template<int D> void testAddition();

SCENARIO("Adding MW trees", "[mw_adder], [tree_builder]") {
    GIVEN("Two MW functions in 1D") {
        testAddition<1>();
    }
    GIVEN("Two MW functions in 2D") {
        testAddition<2>();
    }
    GIVEN("Two MW functions in 3D") {
        testAddition<3>();
    }
}

template<int D> void testAddition() {
    double alpha = 1.0;
    double beta_a = 110.0;
    double beta_b = 50.0;
    double pos_a[3] = {-0.25, 0.35, 1.05};
    double pos_b[3] = {-0.20, 0.50, 1.05};

    GaussFunc<D> a_func(beta_a, alpha, pos_a);
    GaussFunc<D> b_func(beta_b, alpha, pos_b);

    MultiResolutionAnalysis<D> *mra = 0;
    initialize(&mra);

    // Setting up adaptor and TreeBuilders
    double prec = 1.0e-4;
    WaveletAdaptor<D> w_adaptor(prec);
    GridGenerator<D> G(*mra);
    MWProjector<D> Q(*mra, w_adaptor);
    MWAdder<D> add(*mra);

    // Build initial empty grid
    FunctionTree<D> *a_tree = G(a_func);
    FunctionTree<D> *b_tree = G(b_func);

    // Project functions
    Q(*a_tree, a_func);
    Q(*b_tree, b_func);

    // Reference integrals
    const double a_int = a_tree->integrate();
    const double b_int = b_tree->integrate();

    const double a_coef = 1.0;
    const double b_coef = 2.0;

    AdditionVector<D> sum_vec;
    WHEN("the functions are added") {
        sum_vec.push_back(*a_tree, a_coef);
        sum_vec.push_back(*b_tree, b_coef);
        FunctionTree<D> *c_tree = add(sum_vec);
        sum_vec.clear();

        THEN("their integrals add up") {
            double c_int = c_tree->integrate();
            double ref_int = a_coef*a_int + b_coef*b_int;
            REQUIRE( c_int == Approx(ref_int) );
        }

        AND_WHEN("the first function is subtracted") {
            sum_vec.push_back(*c_tree);
            sum_vec.push_back(*a_tree, -1.0);
            FunctionTree<D> *d_tree = add(sum_vec);
            sum_vec.clear();

            THEN("the integral is the same as the second function") {
                double d_int = d_tree->integrate();
                double ref_int = b_coef*b_int;
                REQUIRE( d_int == Approx(ref_int) );
            }

            AND_WHEN("the second function is subtracted") {
                sum_vec.push_back(*d_tree);
                sum_vec.push_back(*b_tree, -b_coef);
                FunctionTree<D> *e_tree = add(sum_vec);
                sum_vec.clear();

                THEN("the integral is zero") {
                    double e_int = e_tree->integrate();
                    double ref_int = 0.0;
                    REQUIRE( e_int == Approx(ref_int) );
                }
                delete e_tree;
            }
            delete d_tree;
        }
        delete c_tree;
    }
    delete b_tree;
    delete a_tree;

    finalize(&mra);
}



} // namespace

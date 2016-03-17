#include "catch.hpp"

#include "factory_functions.h"
#include "GridGenerator.h"
#include "MWProjector.h"
#include "MWMultiplier.h"
#include "WaveletAdaptor.h"

namespace mw_multiplier {

template<int D> void testMultiplication();

SCENARIO("MWMultiplier", "[mw_multiplier], [tree_builder]") {
    GIVEN("Two MW functions in 1D") {
        testMultiplication<1>();
    }
    GIVEN("Two MW functions in 2D") {
        testMultiplication<2>();
    }
    GIVEN("Two MW functions in 3D") {
        testMultiplication<3>();
    }
}

template<int D> void testMultiplication() {
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
    MWMultiplier<D> mult(*mra);

    // Build initial empty grid
    FunctionTree<D> *a_tree = G(a_func);
    FunctionTree<D> *b_tree = G(b_func);

    // Project functions
    Q(*a_tree, a_func);
    Q(*b_tree, b_func);

    MultiplicationVector<D> prod_vec;
    WHEN("the functions are multiplied") {
        prod_vec.push_back(a_tree);
        prod_vec.push_back(b_tree);
        FunctionTree<D> *c_tree = mult(prod_vec);
        prod_vec.clear();
        delete c_tree;
    }
    delete b_tree;
    delete a_tree;

    finalize(&mra);
}



} // namespace

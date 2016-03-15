#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "MWAdder.h"
#include "WaveletAdaptor.h"
#include "CopyAdaptor.h"

namespace mw_adder {

TEST_CASE("MWAdder", "[mw_adder], [tree_builder]") {
    GaussFunc<3> *f_func = 0;
    initialize(&f_func);
    GaussFunc<3> *g_func = 0;
    initialize(&g_func);
	double g_pos[3] = {2.0, 0.5, 0.11};
	g_func->setPos(g_pos);
    MultiResolutionAnalysis<3> *mra = 0;
    initialize(&mra);

    // Setting up adaptor and projector
    double prec = 1.0e-3;
    WaveletAdaptor<3> w_adaptor(prec);
    MWProjector<3> Q(*mra, w_adaptor);
	MWAdder<3> add(*mra, w_adaptor);

	double f_coef = 1.0, g_coef = 2.0;
    FunctionTree<3> *f_tree = Q(*f_func);
    FunctionTree<3> *g_tree = Q(*g_func);
    FunctionTree<3> *h_tree = add(f_coef, *f_tree, g_coef, *g_tree);

    double f_int = f_tree->integrate();
    double g_int = g_tree->integrate();
    double h_int = h_tree->integrate();

	REQUIRE( h_int == Approx(f_coef * f_int + g_coef * g_int) );
}



} // namespace

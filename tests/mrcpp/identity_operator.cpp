#include "catch.hpp"

#include "factory_functions.h"
#include "MWProjector.h"
#include "IdentityKernel.h"

namespace identity_operator {

TEST_CASE("Initialize identity kernel", "[identity_kernel], [identity_operator], [mw_operator]") {
    double prec = 1.0e-4;
    MultiResolutionAnalysis<1> *mra = 0;
    initializeKernel(&mra);
    MWProjector<1> Q(*mra, prec);
    finalize(&mra);

    IdentityKernel id_kern(prec);
    REQUIRE( id_kern.size() == 1 );

    FunctionTree<1> *kern_tree = Q(id_kern);
    REQUIRE( kern_tree->integrate() == Approx(1.0).epsilon(prec) );

    delete kern_tree;
}

} // namespace

#include "catch.hpp"

#include "factory_functions.h"
#include "GreensKernel.h"
#include "PoissonKernel.h"

using namespace std;

namespace poisson_kernel {

TEST_CASE("Poisson's kernel", "[poisson_kernel], [greens_kernel]") {
	
    MultiResolutionAnalysis<1> *mra = 0;
    initialize(&mra);
	PoissonKernel P(*mra);

    SECTION("Default kernel construction") {
		REQUIRE( P.getNTerms() == 44 );
    }

    SECTION("Check some 1/r values") {
		double x = 0.001;
		for (int i = 0; i < 10; i++) {
			REQUIRE( P.evalf(&x) == Approx(1.0/x) );
			x *= 2.0;
		}
    }
}

} // namespace

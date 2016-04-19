#include "catch.hpp"

#include "factory_functions.h"
#include "GreensKernel.h"
#include "HelmholtzKernel.h"

using namespace std;

namespace poisson_kernel {

TEST_CASE("Helmholtz kernel", "[helmholtz_kernel], [greens_kernel]") {
	
    MultiResolutionAnalysis<1> *mra = 0;
    initialize(&mra);
	double mu = 0.1;
	HelmholtzKernel H(mu, *mra);

    SECTION("Default kernel construction") {
		REQUIRE( H.getNTerms() == 46 );
    }

    SECTION("Check some exp(-mu*r)/r values") {
		double x = 0.001;
		for (int i = 0; i < 10; i++) {
			REQUIRE( H.evalf(&x) == Approx(exp(-mu*x)/x) );
			x *= 2.0;
		}
    }
}

} // namespace

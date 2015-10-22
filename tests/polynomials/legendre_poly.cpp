#include "catch.hpp"

#include "LegendrePoly.h"

using namespace std;

SCENARIO("Legendre polynomials can be constructed", "[legendre_constructor], [polynomials]") {
    GIVEN("the the unscaled Legendre polynomial L_5") {
        LegendrePoly L_5(5);
        THEN("the boundaries of L_5 [-1.0, 1.0]") {
            REQUIRE( L_5.getScaledLowerBound() == Approx(-1.0) );
            REQUIRE( L_5.getScaledUpperBound() == Approx(1.0) );
        }
        THEN("the order of L_5 is 5") {
            REQUIRE( L_5.getOrder() == 5 );
        }
        THEN("the function value L_5(1.0) = 1.0") {
            REQUIRE( L_5.evalf(1.0) == Approx(1.0) );
        }
    }
    GIVEN("the Legendre polynomial L_3 shifted to the unit interval") {
        LegendrePoly L_3(3, 2.0, 1.0);
        THEN("the boundaries of L_3 [0.0, 1.0]") {
            REQUIRE( L_3.getScaledLowerBound() == Approx(0.0) );
            REQUIRE( L_3.getScaledUpperBound() == Approx(1.0) );
        }
        THEN("the order of L_k is 3") {
            REQUIRE( L_3.getOrder() == 3 );
        }
        THEN("the function value L_3(1.0) = 1.0") {
            REQUIRE( L_3.evalf(1.0) == Approx(1.0) );
        }
    }
}

SCENARIO("Legendre polynomials are orthogonal", "[legendre_orthogonal], [polynomials]") {
    GIVEN("the set of 10 lowest order Legendre polynomials") {
        vector<LegendrePoly *> L;
        for (int k = 0; k < 10; k++) {
            LegendrePoly *L_k = new LegendrePoly(k, 2.0, 1.0);
            L.push_back(L_k);
        }
        THEN("the overlap <L_0|L_0> is non-zero") {
            LegendrePoly &L_0 = *L[0];
            double overlap = L_0.innerProduct(L_0);
            REQUIRE( fabs(overlap) > MachineZero );
        }
        THEN("the overlap <L_5|L_5> is non-zero") {
            LegendrePoly &L_5 = *L[5];
            double overlap = L_5.innerProduct(L_5);
            REQUIRE( fabs(overlap) > MachineZero );
        }
        THEN("the overlap <L_3|L_2> is zero") {
            LegendrePoly &L_2 = *L[2];
            LegendrePoly &L_3 = *L[3];
            double overlap = L_3.innerProduct(L_2);
            REQUIRE( fabs(overlap) == Approx(0.0) );
        }
        THEN("the overlap <L_1|L_9> is zero") {
            LegendrePoly &L_1 = *L[1];
            LegendrePoly &L_9 = *L[9];
            double overlap = L_1.innerProduct(L_9);
            REQUIRE( fabs(overlap) == Approx(0.0) );
        }
    }
}

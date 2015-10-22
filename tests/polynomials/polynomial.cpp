#include "catch.hpp"

#include <cmath>
#include <iostream>

#include "Polynomial.h"

using namespace Eigen;

SCENARIO("Polynomials can be constructed", "[poly_constructor], [polynomials]") {
    GIVEN("a default polynomial P") {
        Polynomial P;
        THEN("P is the zero function") {
            REQUIRE( P.getOrder() == 0 );
            REQUIRE( P.getCoefs()[0] == Approx(0.0) );
        }
        THEN("P is not dilated") {
            REQUIRE( P.getDilation() == Approx(1.0) );
        }
        THEN("P is not translated") {
            REQUIRE( P.getTranslation() == Approx(0.0) );
        }
        THEN("P is not bounded") {
            REQUIRE_FALSE( P.isBounded() );
        }
    }
    GIVEN("A coefficient vector with trailing zeros") {
        Vector4d c = {0.0, 1.0, 0.0, 0.0};
        WHEN("A polynomial P is constructed") {
            Polynomial P(c);
            THEN("The order of P is given by the highest non-zero coef") {
                REQUIRE( P.getOrder() == 1 );
            }
        }
    }
    GIVEN("The bounded polynomial P(x) = x") {
        double a = 0.0;
        double b = 2.0;
        Vector2d c = {0.0, 1.0};
        Polynomial P(c, &a, &b);
        WHEN("the copy constructor is used, Q(P)") {
            Polynomial Q(P);
            THEN("Q(x) = x, and Q is bounded") {
                REQUIRE( Q.getOrder() == 1 );
                REQUIRE( Q.getCoefs()[0] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(1.0) );
                REQUIRE( Q.getLowerBound(0) == 0.0 );
                REQUIRE( Q.getUpperBound(0) == 2.0 );
            }
        }
        WHEN("the assignment operator is used, Q = P") {
            Polynomial Q;
            Q = P;
            THEN("Q(x) = x, but Q is unbounded") {
                REQUIRE( Q.getOrder() == 1 );
                REQUIRE( Q.getCoefs()[0] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(1.0) );
                REQUIRE_FALSE( Q.isBounded() );
            }
        }
    }
}

SCENARIO("Polynomials can be evaluated", "[poly_evalf], [polynomials]") {
    GIVEN("a bounded polynomial P") {
        double a = 0.0;
        double b = 2.0;
        Vector3d c = {1.0, 0.0, 2.0};
        Polynomial P(c, &a, &b);
        WHEN("P is evaluated within its bounds") {
            double x = 1.5;
            double calc_val = P.evalf(&x);
            THEN("the function value is known") {
                double ref_val = 1.0 + 2.0*x*x;
                REQUIRE( calc_val == Approx(ref_val) );
            }
        }
        WHEN("P is evaluated out of bounds") {
            double x = 2.5;
            double calc_val = P.evalf(&x);
            THEN("the function value is zero") {
                double ref_val = 0.0;
                REQUIRE( calc_val == Approx(ref_val) );
            }
        }
    }
}

SCENARIO("Polynomials can be scaled and translated", "[poly_scale], [polynomials]") {
    GIVEN("a bounded polynomial P") {
        double a = -1.0;
        double b = 1.0;
        Vector3d c = {0.0, 1.0, 1.0};
        Polynomial P(c, &a, &b);
        WHEN("P is rescaled") {
            double n = 2.0;
            double l = 1.0;
            P.rescale(n, l);
            THEN("the dilation changes") {
                REQUIRE( P.getDilation() == Approx(2.0) );
            }
            THEN("the translation changes") {
                REQUIRE( P.getTranslation() == Approx(1.0) );
            }
            THEN("the scaled bounds change") {
                REQUIRE( P.getScaledLowerBound() == Approx(0.0) );
                REQUIRE( P.getScaledUpperBound() == Approx(1.0) );
            }
            THEN("the unscaled bounds don't change") {
                REQUIRE( P.getLowerBound(0) == Approx(-1.0) );
                REQUIRE( P.getUpperBound(0) == Approx(1.0) );
            }
            THEN("the scaled evaluation is known") {
                double x = 0.3;
                double calc_val = P.evalf(&x);
                double ref_val = (2.0*x - 1.0) + (2.0*x - 1.0)*(2.0*x - 1.0);
                REQUIRE( calc_val == Approx(ref_val) );
            }
        }
    }
}

SCENARIO("Polynomials can be added and multiplied", "[poly_arithmetics], [polynomials]") {
    GIVEN("two polynomials P and Q") {
        Vector4d c1 = {0.0, 1.0, 1.0, 0.0};
        Vector2d c2 = {0.0, 1.0};
        Polynomial P(c1);
        Polynomial Q(c2);
        WHEN("Q += P") {
            Q += P;
            THEN("the coefficients of P are unchanged") {
                REQUIRE( P.getOrder() == 2 );
                REQUIRE( P.getCoefs()[0] == Approx(0.0) );
                REQUIRE( P.getCoefs()[1] == Approx(1.0) );
                REQUIRE( P.getCoefs()[2] == Approx(1.0) );
            }
            THEN("the coefficients of Q change") {
                REQUIRE( Q.getOrder() == 2 );
                REQUIRE( Q.getCoefs()[0] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(2.0) );
                REQUIRE( Q.getCoefs()[2] == Approx(1.0) );
            }
        }
        WHEN("R = P + Q") {
            Polynomial R;
            R = P+Q;
            THEN("the coefficients of R are known") {
                REQUIRE( R.getOrder() == 2 );
                REQUIRE( R.getCoefs()[0] == Approx(0.0) );
                REQUIRE( R.getCoefs()[1] == Approx(2.0) );
                REQUIRE( R.getCoefs()[2] == Approx(1.0) );
            }
        }
        WHEN("P *= 2.0") {
            P *= 2.0;
            THEN("the coefficients of P change") {
                REQUIRE( P.getOrder() == 2 );
                REQUIRE( P.getCoefs()[0] == Approx(0.0) );
                REQUIRE( P.getCoefs()[1] == Approx(2.0) );
                REQUIRE( P.getCoefs()[2] == Approx(2.0) );
            }
        }
        WHEN("R = P*3.0") {
            Polynomial R;
            R = P*3.0;
            THEN("the coefficients of R are known") {
                REQUIRE( R.getOrder() == 2 );
                REQUIRE( R.getCoefs()[0] == Approx(0.0) );
                REQUIRE( R.getCoefs()[1] == Approx(3.0) );
                REQUIRE( R.getCoefs()[2] == Approx(3.0) );
            }
        }
        WHEN("Q *= P") {
            Q *= P;
            THEN("the coefficients of P are unchanged") {
                REQUIRE( P.getOrder() == 2 );
                REQUIRE( P.getCoefs()[0] == Approx(0.0) );
                REQUIRE( P.getCoefs()[1] == Approx(1.0) );
                REQUIRE( P.getCoefs()[2] == Approx(1.0) );
            }
            THEN("the coefficients of Q change") {
                REQUIRE( Q.getOrder() == 3 );
                REQUIRE( Q.getCoefs()[0] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[2] == Approx(1.0) );
                REQUIRE( Q.getCoefs()[3] == Approx(1.0) );
            }
        }
        WHEN("R = P * Q") {
            Polynomial R;
            R = P*Q;
            THEN("the coefficients of R are known") {
                VectorXd &cr = R.getCoefs();
                REQUIRE( R.getCoefs()[0] == Approx(0.0) );
                REQUIRE( R.getCoefs()[1] == Approx(0.0) );
                REQUIRE( R.getCoefs()[2] == Approx(1.0) );
                REQUIRE( R.getCoefs()[3] == Approx(1.0) );
            }
        }
    }
}

SCENARIO("Polynomials can be differentiated and integrated", "[poly_diff], [polynomials]") {
    GIVEN("a polynomial P") {
        Vector3d c = {0.0, 1.0, 2.0};
        Polynomial P(c);
        WHEN("Q is the derivative of P") {
            Polynomial Q = P.calcDerivative();
            THEN("the coefficients of Q are known") {
                REQUIRE( Q.getOrder() == 1 );
                REQUIRE( Q.getCoefs()[0] == Approx(1.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(4.0) );
            }
        }
        WHEN("Q is the anti derivative of P") {
            Polynomial Q = P.calcAntiDerivative();
            THEN("the coefficients of Q are known") {
                REQUIRE( Q.getOrder() == 3 );
                REQUIRE( Q.getCoefs()[0] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[1] == Approx(0.0) );
                REQUIRE( Q.getCoefs()[2] == Approx(1.0/2.0) );
                REQUIRE( Q.getCoefs()[3] == Approx(2.0/3.0) );
            }
        }
    }
    GIVEN("The bounded polynomial P on [-1.0, 1.0]") {
        double a = -1.0;
        double b = 1.0;
        Vector3d c = {0.0, 1.0, 1.0};
        Polynomial P(c, &a, &b);
        THEN("P can be integrated on its full domain") {
            double calc_int = P.integrate();
            double ref_int = 2.0/3.0;
            REQUIRE( calc_int == Approx(ref_int) );
        }
        THEN("P can be integrated on a subdomain [0.0, 1.0]") {
            a = 0.0;
            double calc_int = P.integrate(&a, &b);
            double ref_int = 5.0/6.0;
            REQUIRE( calc_int == Approx(ref_int) );
        }
    }
}

SCENARIO("Bounded polynomials have inner products and norms", "[poly_norm], [polynomials]") {
    GIVEN("two unbounded polynomials P and Q") {
        Vector3d c1 = {0.0, 1.0, 1.0};
        Vector2d c2 = {0.0, 1.0};
        Polynomial P(c1);
        Polynomial Q(c2);
        THEN("the norm of P is undefined") {
            REQUIRE( P.calcSquareNorm() < 0.0 );
        }
        WHEN("P is bounded") {
            double a = -1.0;
            double b = 1.0;
            P.setBounds(&a, &b);
            THEN("the inner product <P|Q> is defined") {
                double ref_inner = 2.0/3.0;
                double calc_inner = P.innerProduct(Q);
                REQUIRE( calc_inner == Approx(ref_inner) );
            }
            THEN("the norm of P is defined") {
                double ref_norm = 16.0/15.0;
                double calc_norm = P.calcSquareNorm();
                REQUIRE( calc_norm == Approx(ref_norm) );
            }
            THEN("P can be normalized") {
                P.normalize();
                double calc_norm = P.calcSquareNorm();
                REQUIRE( calc_norm == Approx(1.0) );
            }
        }
    }
}


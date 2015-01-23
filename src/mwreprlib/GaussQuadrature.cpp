/*
 */

#include "GaussQuadrature.h"
#include "LegendrePoly.h"
#include "macros.h"
#include "constants.h"

using namespace Eigen;
using namespace std;

/** Constructor for Gauss-Legendre quadrature.
 *
 * \param order Polynominal order
 * \param a Lower bound of validity
 * \param b Upper bound of validity
 * \param intervals Number of intervals to divde |a-b| into
 */
GaussQuadrature::GaussQuadrature(int order, double a, double b, int intervals) :
    order(order), A(a), B(b), intervals(intervals) {
    if (order < 0 || order > MaxGaussOrder) {
        MSG_ERROR("Gauss quadrature order " << order <<
                " is larger than the maximum of " << MaxGaussOrder);
    }
    if (a >= b) {
        MSG_ERROR("Invalid Gauss interval, a > b.");
    }
    if (intervals < 1) {
        MSG_ERROR("Invalid number of intervals, intervals < 1");
    }
    npts = order * intervals;
    roots = VectorXd::Zero(npts);
    weights = VectorXd::Zero(npts);
    unscaledRoots = VectorXd::Zero(order);
    unscaledWeights = VectorXd::Zero(order);
    // set up unscaled Gauss points and weights ( interval ]-1,1[)
    if (calcGaussPtsWgts() != 1) {
        MSG_ERROR("Setup of Gauss-Legendre weights failed.")
    }
    calcScaledPtsWgts();
}

GaussQuadrature::~GaussQuadrature() {
}

void GaussQuadrature::setBounds(double a, double b) {
    if (fabs(A - a) < MachineZero and fabs(B - b) < MachineZero) {
        return;
    }
    if (a >= b) {
        MSG_ERROR("Invalid bounds: a > b");
    }
    A = a;
    B = b;
    calcScaledPtsWgts();

}

void GaussQuadrature::setIntervals(int i) {
    if (i == intervals) {
        return;
    }
    if (i < 1) {
        MSG_ERROR("Invalid number of integration intervals: " << i);
    }
    intervals = i;
    npts = order * intervals;
    roots = VectorXd::Zero(npts);
    weights = VectorXd::Zero(npts);
    calcScaledPtsWgts();

}

const VectorXd GaussQuadrature::getRoots(double a, double b, int intervals) const {
    if (fabs(A - a) < MachineZero and fabs(B - b) < MachineZero) {
        return roots;
    }
    VectorXd newRoots(order * intervals);
    rescaleRoots(newRoots, a, b, intervals);
    return newRoots;
}

const VectorXd GaussQuadrature::getWeights(double a, double b, int intervals) const {
    if (fabs(A - a) < MachineZero and fabs(B - b) < MachineZero) {
        return weights;
    }
    VectorXd newWeights(order * intervals);
    rescaleWeights(newWeights, a, b, intervals);
    return newWeights;
}

/** Calculate scaled distribution of roots for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::rescaleRoots(VectorXd &rts, double a, double b,
        int intervals) const {

    // lenght of one block
    double transl = (b - a) / (double) intervals;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < order; j++) {
            rts(k) = unscaledRoots(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calculate scaled distribution of weights for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::rescaleWeights(VectorXd &wgts, double a, double b,
        int intervals) const {

    // lenght of one block
    double transl = (b - a) / (double) intervals;

    int k = 0;
    double pos = a;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < order; j++) {
            wgts(k) = unscaledWeights(j) * xl + pos + xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calculate scaled distribution of points and weights for Gauss-Legendre
 * quadrature on on ]a,b[. The number of quadrature points on the interval
 * is scale*(order+1).
 */
void GaussQuadrature::calcScaledPtsWgts() {

    // lenght of one block
    double transl = (B - A) / (double) intervals;

    int k = 0;
    double pos = A;
    double xl = transl * 0.5;
    // scale and translate Gauss points and weights
    for (int i = 0; i < intervals; i++) {
        for (int j = 0; j < order; j++) {
            roots(k) = unscaledRoots(j) * xl + pos + xl;
            weights(k) = unscaledWeights(j) * xl;
            ++k;
        }
        pos = pos + transl;
    }
}

/** Calulate distribution of points and weights for Guass-Legendre quadrature on
 * ]-1,1[.
 *
 * Find quadrature points and weights by solving for the roots of
 * Legendre polynomials using Newtons method. Using double precison the
 * maximum stable order is currently set to 13. Return 1 on success, 0 on failure.
 *
 */
int GaussQuadrature::calcGaussPtsWgts() {
    double z, z1, xm, xl;

    int K;
    if (order % 2 == 0) {
        K = order / 2;
    } else {
        K = (order + 1) / 2;
    }

    double a = -1.0;
    double b = 1.0;

    xm = (b + a) * 0.5;
    xl = (b - a) * 0.5;

    LegendrePoly legendrep(order, 1.0, 0.0); // Interval [-1,1]
    Vector2d lp;

    for (int i = 0; i < K; i++) {
        z = cos(pi * (i + 0.75) / (order + 0.5));
        int iter;
        for (iter = 0; iter < NewtonMaxIter; iter++) {
            lp = legendrep.firstDerivative(z);

            z1 = z;
            z = z1 - lp(0) / lp(1);
            if (fabs(z - z1) <= EPS) {
                break;
            }
        }
        if (iter == NewtonMaxIter) {
            return 0;
        }

        unscaledRoots(i) = xm - xl * z;
        unscaledRoots(order - 1 - i) = xm + xl * z;

        unscaledWeights(i) = 2.e0 * xl / ((1.e0 - z * z) * lp(1) * lp(1));
        unscaledWeights(order - 1 - i) = unscaledWeights(i);
    }
    return 1;
}

/** Integrate a 1D-function f(x) using quadrature */
double GaussQuadrature::integrate(const RepresentableFunction<1> &func) const {

    double isum = 0.e0;
    double r[1];
    for (int i = 0; i < npts; i++) {
        r[0] = roots(i);
        isum += weights(i) * func.evalf(r);
    }
    return isum;
}

/** Integrate a 2D-function f(x1, x2) using quadrature */
double GaussQuadrature::integrate(const RepresentableFunction<2> &func) const {

    double jsum;
    double r[2];
    double isum = 0.e0;
    for (int i = 0; i < npts; i++) {
        jsum = 0.e0;
        r[0] = roots(i);
        for (int j = 0; j < npts; j++) {
            r[1] = roots(j);
            jsum += weights(j) * func.evalf(r);
        }
        isum += jsum * weights(i);

    }
    return isum;
}

/** Integrate a 3D-function f(x1, x2, x3) using quadrature */
double GaussQuadrature::integrate(const RepresentableFunction<3> &func) const {

    double isum, jsum, ksum;
    double r[3];

    isum = 0.e0;
    for (int i = 0; i < npts; i++) {
        jsum = 0.e0;
        r[0] = roots(i);
        for (int j = 0; j < npts; j++) {
            ksum = 0.e0;
            r[1] = roots(j);
            for (int k = 0; k < npts; k++) {
                r[2] = roots(k);
                ksum += weights(k) * func.evalf(r);
            }
            jsum += ksum * weights(j);
        }
        isum += jsum * weights(i);
    }
    return isum;
}

/** Integrate a ND-function f(x1,...), allowing for different
 * quadrature in each dimension.
 *
 * This function has been implemented using a recursive algorithm.
 */
double GaussQuadrature::integrate_nd(const RepresentableFunction<3> &func, int axis) const {
    NEEDS_TESTING

    double sum;
    static double r[MaxQuadratureDim];

    sum = 0.e0;
    for (int i = 0; i < npts; i++) {
        r[axis] = roots(i);
        if (axis < 2) {
            sum += integrate_nd(func, axis) * weights(i);
        } else {
            sum += weights(i) * func.evalf(r);
        }
    }
    return sum;
}


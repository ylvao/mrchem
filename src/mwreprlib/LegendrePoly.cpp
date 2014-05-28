/*
 *
 *  \date Jul 5, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "LegendrePoly.h"
#include "constants.h"
#include "macros.h"

using namespace std;
using namespace Eigen;

typedef ObjectCache<LegendrePoly> LegendreCache;

/** Legendre polynomial constructed on [-1,1] and scaled by n and translated by l */
LegendrePoly::LegendrePoly(int _order, double n, double l) :
    Polynomial(_order) {
    // Since we create Legendre polynomials recursively on [-1,1]
    // we cache all lower order polynomilas for future use.
    LegendreCache &Cache = LegendreCache::getInstance();
    if (_order >= 1) {
	if (not Cache.hasId(_order - 1)) {
	    LegendrePoly *lp = new LegendrePoly(_order - 1);
	    Cache.load(_order - 1, lp, 2 * sizeof(double) * (_order + 1));
	}
    }
    computeLegendrePolynomial(_order);
    haveWeightsAndRoots = false;
    double a = -1.0;
    double b = 1.0;
    setBounds(&a, &b);
    rescale(n, l);
}

/** Compute Legendre polynomial coefs on interval [-1,1] */
void LegendrePoly::computeLegendrePolynomial(int _order) {
    assert(this->size() >= _order);
    if (_order == 0) {
	this->coefs(0) = 1.0;
    } else if (_order == 1) {
	this->coefs(0) = 0.0;
	this->coefs(1) = 1.0;
    } else {
	LegendreCache &Cache = LegendreCache::getInstance();
	LegendrePoly &p1 = Cache.get(_order - 1);
	LegendrePoly &p2 = Cache.get(_order - 2);

	this->coefs(0) = -((1.0 * _order - 1.0) * p2.coefs(0)) / (1.0 * _order);
	for (int j = 1; j < _order + 1; j++) {
	    if (j <= _order - 2) {
		this->coefs(j) = (((2.0 * _order - 1.0) * p1.coefs(j - 1)) - ((1.0
				* _order - 1.0) * p2.coefs(j))) / (1.0 * _order);
	    } else {
		this->coefs(j) = (((2.0 * _order - 1.0) * p1.coefs(j - 1))) / (1.0 * _order);
	    }
	}
    }
}

/** Calculate the value of an n:th order Legendre polynomial in x.
 * This routine is not used, except for debugging, since it does not
 * conform to the general polynomial class structure.
 *
 * @param x Pretty obvious, right?
 * @return Value of Pn(x)
 */
double LegendrePoly::evalLegendrePoly(double x) const {
    NOT_IMPLEMENTED_ABORT
	double c1, c2, c4, ym, yp, y;

	double q = this->N * x + this->L;
	if (out_of_bounds(x, this->A[0], this->B[0])) {
		MSG_FATAL("Argument out of bounds: " << x << " [" <<
				this->A[0] << ", " << this->B[0] << "]");
	}
	int order = getOrder();
	if (order == 0)
		return 1.e0;
	if (order == 1)
		return q;

	y = q;
	yp = 1.e0;
	for (int i = 2; i < order + 1; i++) {
		c1 = (double) i;
		c2 = c1 * 2.e0 - 1.e0;
		c4 = c1 - 1.e0;
		ym = y;
		y = (c2 * q * y - c4 * yp) / c1;
		yp = ym;
	}
	return y;
}

/** Calulate the value of an n:th order Legendre polynominal in x, including
 * the first derivative.
 */

Vector2d LegendrePoly::firstDerivative(double x) const {
    double c1, c2, c4, ym, yp, y;
    double dy, dyp, dym;

    if (out_of_bounds(x, this->A[0], this->B[0])) {
	MSG_FATAL("Argument out of bounds: " << x << " [" <<
		this->A[0] << ", " << this->B[0] << "]");
    }

    double q = this->N * x + this->L;
    Vector2d val;

    int order = getOrder();
    if (order == 0) {
	val(0) = 1.0;
	val(1) = 0.0;
	return val;
    }

    if (order == 1) {
	val(0) = q;
	val(1) = this->N * 1.0 + this->L;
	return val;
    }

    y = q;
    dy = 1.0;
    yp = 1.0;
    dyp = 0.0;
    for (int i = 2; i < order + 1; i++) {
	c1 = (double) i;
	c2 = c1 * 2.0 - 1.0;
	c4 = c1 - 1.0;
	ym = y;
	y = (c2 * q * y - c4 * yp) / c1;
	yp = ym;
	dym = dy;
	dy = (c2 * q * dy - c4 * dyp + c2 * yp) / c1;
	dyp = dym;
    }

    val(0) = y;
    val(1) = dy;
    return val;
}

/** Calulate the value of an n:th order Legendre polynominal in x, including
 * first and second derivatives.
 */
Vector3d LegendrePoly::secondDerivative(double x) const {
    NOT_IMPLEMENTED_ABORT
	double c1, c2, c4, ym, yp, y, d2y;
	double dy, dyp, dym, d2ym, d2yp;

	double q = this->N * x + this->L;
	if (out_of_bounds(x, this->A[0], this->B[0])) {
		MSG_FATAL("Argument out of bounds: " << x << " [" <<
				this->A[0] << ", " << this->B[0] << "]");
	}

	Vector3d val;

	int order = getOrder();
	if (order == 0) {
		val(0) = 1.e0;
		val(1) = 0.e0;
		val(2) = 0.e0;
		return val;
	}

	if (order == 1) {
		val(0) = q;
		val(1) = this->N * 1.e0 + this->L;
		val(2) = 0.e0;
		return val;
	}

	y = q;
	dy = 1.e0;
	d2y = 0.e0;
	yp = 1.e0;
	dyp = 0.e0;
	d2yp = 0.e0;
	for (int i = 2; i < order + 1; i++) {
		c1 = (double) i;
		c2 = c1 * 2.e0 - 1.e0;
		c4 = c1 - 1.e0;
		ym = y;
		y = (c2 * x * y - c4 * yp) / c1;
		yp = ym;
		dym = dy;
		dy = (c2 * x * dy - c4 * dyp + c2 * yp) / c1;
		dyp = dym;
		d2ym = d2y;
		d2y = (c2 * x * d2y - c4 * d2yp + c2 * 2.e0 * dyp) / c1;
		d2yp = d2ym;
	}
	val(0) = y;
	val(1) = dy;
	val(2) = d2y;
	return val;
}

/** Compute the roots of the Legendre polynomial of
 degree K together with the corresponding weights, on interval [a,b]. The code here is
 a rewrite of the Matlab code in file: lgwt.m.*/
/*N=N-1;
 // N1=N+1; N2=N+2;
 xu=linspace(-1,1,N1)';
 % Initial guess
 y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
 % Legendre-Gauss Vandermonde Matrix
 L=zeros(N1,N2);
 % Derivative of LGVM
 Lp=zeros(N1,N2);
 % Compute the zeros of the N+1 Legendre Polynomial
 % using the recursion relation and the Newton-Raphson method
 y0=2 */
void LegendrePoly::calcRootsAndWeights(double prec) {
    NOT_IMPLEMENTED_ABORT

	//TODO: compare to algo in GaussQuadrature!
	if (prec < 0.0L) {
		MSG_FATAL("Requested precision < 0.0");
	}
	int order = getOrder();
	roots = VectorXd(order);
	weights = VectorXd(order);

	int K = order;
	int N1 = K - 1;
	int K1 = K + 1;
	VectorXd xu(K); //Matrix<double, K, 1>
	VectorXd y(K);
	VectorXd y0(K);
	MatrixXd L(K, K1);
	VectorXd Lp(K);

	L.setZero();

	for (int i = 0; i < K; i++) {
		xu(i) = -1.0 + (2.0 * i) / N1;
		y(i) = cos(((2 * i + 1) * pi) / (2 * N1 + 2)) + ((0.27L) / K) * sin((pi
				* xu(i) * N1) / K1);
	}
	y0.setConstant(2.0);

	double maxDiff = prec + 1.0;
	while (maxDiff > prec) {

		L.col(0).setOnes();
		L.col(1) = y;
		for (int i = 1; i < K; i++) {
			L.col(i + 1) = y.array() * L.col(i).array();
			L.col(i + 1) *= (2.0 * i - 1.0);
			L.col(i + 1) += (1.0 * i - 1.0) * L.col(i - 1);
			L.col(i + 1) *= 1.0 / i;
		}

		Lp = K1 * (L.col(K).array() - y.array() * L.col(K1).array());
		y0 = VectorXd::Ones(K).array() - y.array() * y.array();
		Lp = Lp.array() / y0.array();
		y0 = y;
		y = y0.array() - L.col(K).array() / Lp.array();

		maxDiff = 0.0L;
		for (int i = 0; i < K; i++) {
			maxDiff = max(maxDiff, fabs(y(i) - y0(i)));
		}
	}

	/*Linear map from[-1,1] to [a,b]:
	 x=(a*(1-y)+b*(1+y))/2;
	 Compute the weights:
	 w=(b-a)./((1-y.^2).*Lp.^2)*(K1/K)^2; */
	for (int i = 0; i < K; i++) {
		this->roots(K - 1 - i) = 0.5 * (this->A[0] * (1.0 - y(i)) + this->B[0] * (1.0
				+ y(i)));
		this->weights(K - 1 - i) = (this->B[0] - this->A[0]) * (K1 * K1) / (K * K);
		this->weights(K - 1 - i) /= ((1.0 - y(i) * y(i)) * Lp(i) * Lp(i));
	}
	this->haveWeightsAndRoots = true;
}

/** Change the weights and roots according to the transformation of variable:
 * x->n*x+l. */
void LegendrePoly::rescale(double n, double l) {
    Polynomial::rescale(n, l);
    if (haveWeightsAndRoots) {
	for (int k = 0; k < this->roots.size(); k++) {
	    this->roots(k) = (this->roots(k) - this->L) / this->N;
	    this->weights(k) *= 1.0/this->N;
	}
    }
}


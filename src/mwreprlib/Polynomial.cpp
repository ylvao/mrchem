/**
 *
 * \date Jun 7, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 */

#include <cmath>

#include "macros.h"
#include "Polynomial.h"
#include "MultiException.h"
#include "MathUtils.h"
#include "QuadratureCache.h"
#include "GaussQuadrature.h"


using namespace Eigen;
using namespace std;

Polynomial::Polynomial(int power, const double *a, const double *b) 
	: RepresentableFunction(a, b) {
    assert(power >= 0);
    this->N = 1.0;
    this->L = 0.0;
    this->squareNorm = -1.0;
    this->coefs = VectorXd::Zero(power + 1);
}

/** Makes a complete copy of the polynomail */
Polynomial::Polynomial(const Polynomial &poly) : RepresentableFunction<1>(poly) {
    this->N = poly.getDilation();
    this->L = poly.getTranslation();
    this->coefs = poly.getCoefs();
    this->squareNorm = poly.getSquareNorm();
}

/** Copies only the function, not its bounds */
Polynomial& Polynomial::operator=(const Polynomial &poly) {
    RepresentableFunction<1>::operator=(poly);
    this->N = poly.getDilation();
    this->L = poly.getTranslation();
    this->coefs = poly.getCoefs();
    this->squareNorm = -1.0;
    return *this;
}

double Polynomial::evalf(double x) const {
    if (this->outOfBounds(&x)) {
	return 0.0;
    }
    double xp = 1.0;
    double y = 0.0;
    for (int k = 0; k < getOrder() + 1; k++) {
        y += (xp * this->coefs(k));
        xp *= this->N * x - (double) this->L;
    }
    return y;
}

double Polynomial::evalf(const double *r) const {
    return evalf(*r);
}

/** This returns the actual scaled lower bound */
double Polynomial::getLowerBound(int i) const {
    if (not isBounded()) {
	THROW_ERROR("Unbounded polynomial");
    }
    return (1.0/this->N * (this->A[0] + (double) this->L));
}

/** This returns the actual scaled upper bound */
double Polynomial::getUpperBound(int i) const {
    if (not isBounded()) {
	THROW_ERROR("Unbounded polynomial");
    }
    return (1.0/this->N * (this->B[0] + (double) this->L));
}

/** This returns the original unscaled lower bounds */
const double* Polynomial::getLowerBounds() const {
    return this->A;
}

/** This returns the original unscaled upper bounds */
const double* Polynomial::getUpperBounds() const {
    return this->B;
}

void Polynomial::setCoefs(const VectorXd &c) {
    this->coefs = c;
    this->squareNorm = -1.0;
}

void Polynomial::normalize() {
    if (getSquareNorm() < 0.0) {
	calcSquareNorm();
    }
    double norm = sqrt(getSquareNorm());
    multInPlace(1.0/norm);
    calcSquareNorm();
}

/** Compute the squared L2-norm of Polynomial P(nx+l) over
 * interval [A,B], using Gauss-Legendre quadratue with weight function
 * identical to 1. The L2-norm is stored internally. */
void Polynomial::calcSquareNorm() {
    if (not this->isBounded()) {
        MSG_ERROR("Cannot calculate norm of unbounded polynomials");
    }
    double norm = this->innerProduct(*this);
    if (norm < 0.0) {
        THROW_ERROR("Undefined, L2-norm < 0.0: " << squareNorm);
    }
    if (norm < MachineZero) {
        norm = 0.0;
    }
    this->squareNorm = norm;
}

/** Dilates the polynomial, keeps the domain [A,B] */
void Polynomial::setDilation(double n) {
    if (n <= 0.0) {
	THROW_ERROR("Scaling factor must be positive");
    }
    this->N = n;
    this->squareNorm = -1.0;
}

/** Translates the polynomial, keeps the domain [A,B] */
void Polynomial::setTranslation(double l) {
    this->L = l;
    this->squareNorm = -1.0;
}
/** Dilates and translates the polynomial, changes the domain [A,B] accordingly
  * Transform: P(2^(-n)*x+l)->P(2^(-n')*(2^(-n)x+l)+l') for given arguments n,l. */
void Polynomial::rescale(double n, double l) {
    if (n <= 0.0) {
	THROW_ERROR("Scaling factor must be positive " << n);
    }
    this->L = this->L + l;
    this->N *= n;
    this->squareNorm = -1.0;
}

/** Returns the order of the highest non-zero coef, not the length of the coefs vector */
int Polynomial::getOrder() const {
    int n = 0;
    for (int i = 0; i < this->coefs.size(); i++) {
	if (fabs(this->coefs[i]) > MachineZero) {
	    n = i;
	}
    }
    return n;
}

void Polynomial::clearCoefs() {
    this->coefs = VectorXd::Zero(1);
    this->squareNorm = -1.0;
}

void Polynomial::setZero() {
    int n = this->coefs.size();
    this->coefs = VectorXd::Zero(n);
}

/** Calculate P = c*P */
void Polynomial::multInPlace(double c) {
    Polynomial &P = *this;
    P.setCoefs(c*P.getCoefs());
}

/** Calculate P = P*Q */
void Polynomial::multInPlace(const Polynomial &Q) {
    Polynomial &P = *this;
    if (fabs(P.getDilation() - Q.getDilation()) > MachineZero) {
	THROW_ERROR("Polynomials not defined on same scale.");
    }
    if (fabs(P.getTranslation() - Q.getTranslation()) > MachineZero) {
	THROW_ERROR("Polynomials not defined on same translation.");
    }
    
    int P_order = P.getOrder();
    int Q_order = Q.getOrder();
    int new_order = P_order + Q_order;
    VectorXd coefs = VectorXd::Zero(new_order + 1);
    for (int i = 0; i < P_order + 1; i++) {
        for (int j = 0; j < Q_order + 1; j++) {
            coefs(i + j) += P.coefs(i) * Q.coefs(j);
        }
    }
    P.setCoefs(coefs);
}

/** Calculate Q = c*P */
Polynomial Polynomial::mult(double c) const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.multInPlace(c);
    return Q;
}

/** Calculate R = P*Q. Returns unbounded polynomial. */
Polynomial Polynomial::mult(const Polynomial &Q) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R.multInPlace(Q);
    return R;
}

/** Calculate P = P + c*Q, with a default c = 1.0 */
void Polynomial::addInPlace(const Polynomial &Q, double c) {
    Polynomial &P = *this;
    if (fabs(P.getDilation() - Q.getDilation()) > MachineZero) {
	THROW_ERROR("Polynomials not defined on same scale.");
    }
    if (fabs(P.getTranslation() - Q.getTranslation()) > MachineZero) {
	THROW_ERROR("Polynomials not defined on same translation.");
    }
    
    int P_order = P.getOrder();
    int Q_order = Q.getOrder();
    int new_order = max(P_order, Q_order);
    VectorXd coefs = VectorXd::Zero(new_order + 1);

    for (int i = 0; i < new_order + 1; i++) {
	if (i <= P_order) {
	    coefs[i] += P.getCoefs()[i];
	}
	if (i <= Q_order) {
	    coefs[i] += c*Q.getCoefs()[i];
	}
    }
    P.setCoefs(coefs);
}

/** Calculate R = P + c*Q, with a default c = 1.0. Returns unbounded polynomial. */
Polynomial Polynomial::add(const Polynomial &Q, double c) const {
    const Polynomial &P = *this;
    Polynomial R;
    R = P;
    R.addInPlace(Q, c);
    return R;
}

/** Calculate Q = dP/dx */
Polynomial Polynomial::calcDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcDerivativeInPlace();
    return Q;
}

/** Calculate P = dP/dx */
void Polynomial::calcDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order);
    for (int i = 0; i < newCoefs.size(); i++) {
	newCoefs[i] = double (i+1) * oldCoefs[i+1];
    }
    P.setCoefs(newCoefs);
}

/** Calculate indefinite integral Q = \int dP dx, integration constant set to zero */
Polynomial Polynomial::calcAntiDerivative() const {
    const Polynomial &P = *this;
    Polynomial Q(P);
    Q.calcAntiDerivativeInPlace();
    return Q;
}

/** Calculate indefinite integral P = \int dP dx, integration constant set to zero */
void Polynomial::calcAntiDerivativeInPlace() {
    Polynomial &P = *this;
    int P_order = P.getOrder();
    const VectorXd &oldCoefs = P.getCoefs();
    VectorXd newCoefs = VectorXd::Zero(P_order + 2);
    newCoefs[0] = 0.0;
    newCoefs[1] = oldCoefs[0];
    for (int i = 2; i < newCoefs.size(); i++) {
	newCoefs[i] = 1.0/i * oldCoefs[i-1];
    }
    P.setCoefs(newCoefs);
}

/** Integrate the polynomial P on [a,b] analytically */
double Polynomial::integrate(const double *a, const double *b) const {
    double lb, ub;
    if (a == 0) {
	if (not this->isBounded()) {
	    THROW_ERROR("Cannot integrate polynomial without bounds");
	}
	lb = getLowerBound();
    } else {
	if (this->outOfBounds(a)) {
	    THROW_ERROR("Integration out of bounds");
	}
	lb = a[0];
    }
    if (b == 0) {
	if (not this->isBounded()) {
	    THROW_ERROR("Cannot integrate polynomial without bounds");
	}
	ub = getUpperBound();
    } else {
	if (this->outOfBounds(b)) {
	    THROW_ERROR("Integration out of bounds");
	}
	ub = b[0];
    }
    double sfac = 1.0/this->N;
    Polynomial antidiff = calcAntiDerivative();
    return sfac*(antidiff.evalf(ub) - antidiff.evalf(lb));
}

/* Compute <P,Q> analytically on interval defined by the calling polynomial */
double Polynomial::innerProduct(const Polynomial &Q) const {
    const Polynomial &P = *this;
    if (not P.isBounded()) {
	THROW_ERROR("Unbounded polynomial");
    }
    Polynomial pq = P.mult(Q);
    pq.setBounds(P.getLowerBounds(), P.getUpperBounds());
    return pq.integrate();
}

/** Compute <P,Q> using quadrature. This MIGHT be less accurate than
 the true polynomial product.*/
double Polynomial::innerProduct(const RepresentableFunction<1> &Q, int quadOrder) const {
    const Polynomial &P = *this;
    NOT_IMPLEMENTED_ABORT;
    //FunctionProduct<1, SeparableFunction<1> > fprod(this, &Q);
    //getQuadratureCache(qCache);
    //GaussQuadrature &quad = qCache.get(quadratureOrder);
    //return quad.integrate(fprod);
}
/*
Polynomial Polynomial::projectPolynomial(Polynomial &q) {
    double proj = calcPolynomialProjection(q);
    return (*this) * proj;
}

Polynomial Polynomial::projectFunction(RepresentableFunction<1>&q) {
    double proj = calcFunctionProjection(q);
    return (*this) * proj;
}
*/

/** Compute the projection of P on Q over interval [A,B],
 using Gauss-Legendre quadrature. */
/*
double Polynomial::calcPolynomialProjection(Polynomial &q) {
    double norm_p = getSquareNorm();
    double norm_q = q.getSquareNorm();

    if (norm_p < SQUARE(MachinePrec) || norm_q < SQUARE(MachinePrec)) {
        return 0.0;
    }
    double dot_p = polynomialInnerproduct(q);
    if (norm_p < SQUARE(dot_p) / norm_q) {
        THROW_ERROR("Projection violates the Schwartz inequality; Result is undefined.");
    }
    return dot_p/norm_p;
}
*/
/** Compute the projection of P on a generic Q over interval [A,B],
 using Gauss-Legendre quadrature. This MIGHT be less accurate than
 the true polynomial product. */
/*
double Polynomial::calcFunctionProjection(RepresentableFunction<1> &q) {
    double norm_p = getSquareNorm();
    double dot_p = functionInnerproduct(q);
    return dot_p / norm_p;
}
*/

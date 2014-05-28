/**
 *
 * \date Jun 7, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *  Base class for general polynomials with reasonably advanced
 * properties. The Polynomial class(es) are not implemented in the
 * most efficient manner, because they are only evaluated a fixed
 * number of times in a few predefined points, and all other
 * evaluations are done by linear transformations. PolynomialCache
 * implements the fast, and static const versions of the various
 * 4Polynomials.
 */

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <Eigen/Core>

#include "RepresentableFunction.h"
#include "GaussQuadrature.h"

class Polynomial: public RepresentableFunction<1> {
public:
    Polynomial(int power = 0, const double *a = 0, const double *b = 0);
    Polynomial(const Polynomial &poly);
    Polynomial &operator=(const Polynomial &poly);
    virtual ~Polynomial() {}

    double evalf(double x) const;
    double evalf(const double *r) const;

    double getUpperBound(int i = 0) const; ///< Scaled lower bound
    double getLowerBound(int i = 0) const; ///< Scaled upper bound

    const double *getUpperBounds() const; ///< Unscaled lower bounds
    const double *getLowerBounds() const; ///< Unscaled lower bounds

    void normalize();
    void calcSquareNorm();
    double getSquareNorm() const { return this->squareNorm; }

    double getTranslation() const { return this->L; }
    double getDilation() const { return this->N; }

    void setDilation(double n);
    void setTranslation(double l);
    virtual void rescale(double n, double l);

    int size() const { return this->coefs.size(); } ///< Length of coefs vector
    int getOrder() const; ///< Polynomial order based on highest nonzero coef
    void clearCoefs(); ///< Re-initialize coefs vector to include a single zero entry
    void setZero(); ///< Sets all coefs to zero, keeps length of coefs vector
    void setCoefs(const Eigen::VectorXd &c); ///< Sets the coefs vector and redefines the polynomial order

    Eigen::VectorXd &getCoefs() { return this->coefs; }
    const Eigen::VectorXd &getCoefs() const { return this->coefs; }

    void multInPlace(double c);
    void multInPlace(const Polynomial &p);
    void addInPlace(const Polynomial &p, double c = 1.0);

    Polynomial mult(double c) const;
    Polynomial mult(const Polynomial &p) const;
    Polynomial add(const Polynomial &p, double c = 1.0) const;

    Polynomial calcDerivative() const;
    Polynomial calcAntiDerivative() const;

    void calcDerivativeInPlace();
    void calcAntiDerivativeInPlace();

    double integrate(const double *a = 0, const double *b = 0) const;
    double innerProduct(const Polynomial &p) const;
    double innerProduct(const RepresentableFunction<1> &p, int quadOrder = -1) const;

    //double calcPolynomialProjection(Polynomial &p, int quadOrder = -1) const;
    //double calcFunctionProjection(RepresentableFunction<1> &p, int quadOrder = -1) const;

    //double projectPolynomial(Polynomial &p, int quadOrder = -1);
    //double projectFunction(RepresentableFunction<1> &p, int quadOrder = -1);

    Polynomial operator*(double c) { return mult(c); }
    Polynomial operator*(const Polynomial &p) { return mult(p); }
    Polynomial operator+(const Polynomial &p) { return add(p, 1.0); }
    Polynomial operator-(const Polynomial &p) { return add(p, -1.0); }
    Polynomial &operator*=(double c) { this->multInPlace(c); return *this; }
    Polynomial &operator*=(const Polynomial &p) {multInPlace(p);return *this;}
    Polynomial &operator+=(const Polynomial &p) {addInPlace(p,1.0);return *this;}
    Polynomial &operator-=(const Polynomial &p) {addInPlace(p,-1.0);return *this;}
protected:
    double N; ///< Dilation coeff
    double L; ///< Translation coeff
    double squareNorm;
    Eigen::VectorXd coefs; ///< Expansion coefficients
};
#endif /* POLYNOMIAL_H_ */

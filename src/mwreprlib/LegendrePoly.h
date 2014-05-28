/*
 * 
 *
 *  \date Jul 5, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef LEGENDREPOLY_H_
#define LEGENDREPOLY_H_

#include <Eigen/Core>

#include "Polynomial.h"
#include "ObjectCache.h"

class LegendrePoly: public Polynomial {
public:
	LegendrePoly(int _order, double n=1.0, double l=0.0);
	virtual ~LegendrePoly() {}

	double evalLegendrePoly(double x) const;

	Eigen::Vector2d firstDerivative(double x) const;
	Eigen::Vector3d secondDerivative(double x) const;

	void rescale(double n, double l);
	void calcRootsAndWeights(double prec=1.0e-6);
protected:
	bool haveWeightsAndRoots;
	void computeLegendrePolynomial(int _order);
	Eigen::VectorXd roots;
	Eigen::VectorXd weights;

};

#endif /* LEGENDREPOLY_H_ */

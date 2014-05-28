/*
 * 
 *
 *  \date Jul 22, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef GAUSSQUADRATURE_H_
#define GAUSSQUADRATURE_H_

#include <Eigen/Core>

#include "RepresentableFunction.h"

const int MaxGaussOrder = 42;
static const double EPS = 3.0e-12;
static const int NewtonMaxIter = 10;
static const int MaxQuadratureDim = 7;

class GaussQuadrature {
public:
	GaussQuadrature(int order, double a = -1.0, double b = 1.0,
			int intervals = 1);
	virtual ~GaussQuadrature();
	double integrate(const RepresentableFunction<1> &func) const;
	double integrate(const RepresentableFunction<2> &func) const;
	double integrate(const RepresentableFunction<3> &func) const;
	void setIntervals(int i);
	void setBounds(double a, double b);
	int getIntervals() const {
		return intervals;
	}
	double getUpperBound() const {
		return B;
	}
	double getLowerBound() const {
		return A;
	}
	const Eigen::VectorXd &getRoots() const {
		return roots;
	}
	const Eigen::VectorXd &getWeights() const {
		return weights;
	}
	const Eigen::VectorXd getRoots(double a, double b, int intervals = 1) const;
	const Eigen::VectorXd
			getWeights(double a, double b, int intervals = 1) const;
	const Eigen::VectorXd &getUnscaledRoots() const {
		return unscaledRoots;
	}
	const Eigen::VectorXd &getUnscaledWeights() const {
		return unscaledWeights;
	}
protected:
	int order;
	double A;
	double B;
	int intervals;
	int npts;
	Eigen::VectorXd roots;
	Eigen::VectorXd weights;
	Eigen::VectorXd unscaledRoots;
	Eigen::VectorXd unscaledWeights;
	void rescaleRoots(Eigen::VectorXd &roots, double a, double b,
			int intervals = 1) const;
	void rescaleWeights(Eigen::VectorXd &weights, double a, double b,
			int intervals = 1) const;
	void calcScaledPtsWgts();
	int calcGaussPtsWgts();
	double integrate_nd(const RepresentableFunction<3> &func, int axis = 0) const;
};

#endif /* GAUSSQUADRATURE_H_ */

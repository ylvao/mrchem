/** 
 *
 * \date Jun 6, 2009
 * \author Jonas Juselius <jonas.juselius@uit.no> \n
 *         CTCC, University of Tromsø
 *
 *
 *
 * - Monodimensional gaussian expansion:
 *
 * \f$ g(x) = \sum_{i=1}^M g_i(x_i)
 * = \sum_{i=1}^M \alpha_i e^{-c_i (x-x^0)^2} \f$
 *
 * - Multidimensional gaussian expansion:
 *
 * \f$ G(x) = \sum_{j=1}^M G_j(x)
 * = \sum_{j=1}^M \prod_{i=1}^d g_ij(x_i)
 * = \sum_{j=1}^M \prod_{i=1}^d \alpha_{ij} e^{-c_{ij} (x_i-x_i^0)^2} \f$
 *
 *
 */

#ifndef GAUSSEXP_H_
#define GAUSSEXP_H_

#include <vector>

#include "Gaussian.h"
#include "GaussFunc.h"
#include "NodeIndex.h"

#define GAUSS_EXP_PREC 1.e-10

template<int D>
class GaussExp: public RepresentableFunction<D> {
public:
	GaussExp(int nTerms = 0, double prec = GAUSS_EXP_PREC);
	GaussExp(const GaussExp<D> &gExp);
	GaussExp(const GaussPoly<D> &gPoly);
	GaussExp &operator=(const GaussExp<D> &gExp);
	~GaussExp() { }

	double calcCoulombEnergy();
	double calcSquareNorm();
	void normalize();

	void calcScreening(double nStdDev = defaultScreening);
	bool checkSeedNode(MWNode<D> &node);

	double evalf(const double *r) const;

	GaussExp<D> differentiate(int dir);

	GaussExp<D> add(GaussExp<D> &g);
	GaussExp<D> add(Gaussian<D> &g);
	GaussExp<D> mult(GaussExp<D> &g);
	GaussExp<D> mult(GaussFunc<D> &g);
	GaussExp<D> mult(GaussPoly<D> &g);
	GaussExp<D> mult(double d);
	void multInPlace(double d);

	GaussExp<D> operator+(GaussExp<D> &g) { return this->add(g); }
	GaussExp<D> operator+(Gaussian<D> &g) { return this->add(g); }
	GaussExp<D> operator*(GaussExp<D> &g) { return this->mult(g); }
	GaussExp<D> operator*(GaussFunc<D> &g) { return this->mult(g); }
	GaussExp<D> operator*(GaussPoly<D> &g) { return this->mult(g); }
	GaussExp<D> operator*(double d) { return this->mult(d); }
	void operator*=(double d) { this->multInPlace(d); }

	double getScreening() const { return screening; }
	double getExp(int i) const { return this->funcs[i]->getExp(); }
	double getCoef(int i) const { return this->funcs[i]->getCoef(); }
	const int *getPower(int i) const { return this->funcs[i]->getPower(); }
	const double *getPos(int i) const { return this->funcs[i]->getPos(); }

	double getSquareNorm() {
		if (squareNorm < 0) {
			calcSquareNorm();
		}
		return squareNorm;
	}

	void setGaussFunc(int i, const GaussPoly<D> &g, double c = 1.0);
	void setGaussFunc(int i, const GaussFunc<D> &g, double c = 1.0);

	void setScreen(bool screen);
	void setExp(int i, double a) { this->funcs[i]->setExp(a); }
	void setCoef(int i, double b) { this->funcs[i]->setCoef(b); }
	void setPower(int i, const int power[D]) { this->funcs[i]->setPower(power); }
	void setPos(int i, const double pos[D]) { this->funcs[i]->setPos(pos); }

	void setDefaultScreening(double screen) {
		if (screen < 0) {
			MSG_ERROR("Screening constant cannot be negative!");
		}
		defaultScreening = screen;
	}

	void append(const Gaussian<D> &g);
	void append(const GaussExp<D> &g);

	void calcWaveletCoefs(MWNode<D> &node);

	friend std::ostream& operator<< (std::ostream &o, const GaussExp &gExp) {
		o << "Gaussian Expansion: " << gExp.getNFuncs() << " terms" << std::endl;
		for (int i = 0; i < gExp.getNFuncs(); i++) {
			o << "Term " << i << ":" << std::endl;
			o << gExp.getFunc(i) << std::endl << std::endl;
		}
		return o;
	}
protected:
	static double defaultScreening;
	double screening;
	double squareNorm;
};

#endif /* GAUSSEXP_H_ */

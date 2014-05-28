/*
 * 
 *
 *  \date Sep 27, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 *  Base class of functions that is representable in the mw basis.
 * This includes gaussians, expansions, polynomials and even function trees.
 */

#ifndef REPRESENTABLEFUNCTION_H_
#define REPRESENTABLEFUNCTION_H_

#include "TelePrompter.h"
#include "constants.h"

template<int D>
class RepresentableFunction {
public:
    RepresentableFunction();
    RepresentableFunction(const double *a, const double *b);
    RepresentableFunction(const RepresentableFunction<D> &func);
    RepresentableFunction<D> &operator=(const RepresentableFunction<D> &func);
    virtual ~RepresentableFunction();

    virtual double evalf(const double *r) const = 0;

    void setBounds(const double *a, const double *b);
    void clearBounds();

    bool isBounded() const { return this->bounded; }
    bool outOfBounds(const double *r) const;

    virtual double getUpperBound(int i = 0) const { return this->B[i]; }
    virtual double getLowerBound(int i = 0) const { return this->A[i]; }

    virtual const double *getUpperBounds() const { return this->B; }
    virtual const double *getLowerBounds() const { return this->A; }

    virtual void plotCleanup() { }

    virtual bool isVisibleAtScale(int scale, int nQuadPts) const { return true; }
    virtual bool isZeroOnInterval(const double *a, const double *b) const { return false; }

    friend std::ostream& operator<<(std::ostream &o, const RepresentableFunction<D> &func) {
	o << "RepresentableFunction: " << std::endl;
	o << "  A=[ ";
	for (int i = 0; i < D; i++) {
	    o << func.A[i] << " ";
	}
	o << "]" << std::endl;
	o << "  B=[ ";
	for (int i = 0; i < D; i++) {
	    o << func.B[i] << " ";
	}
	o << "]" << std::endl;
	return o;
    }
protected:
    bool bounded;
    double *A; ///< Lower bound, NULL if unbounded
    double *B; ///< Upper bound, Null if unbounded
};

#endif /* REPRESENTABLEFUNCTION_H_ */

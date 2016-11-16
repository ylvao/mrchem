#ifndef XCPOTENTIAL_H
#define XCPOTENTIAL_H

#include "XCOperator.h"

class XCFunctional;
template<int D> class FunctionTreeVector;

class XCPotential : public XCOperator {
public:
    XCPotential(double prec, XCFunctional &func, OrbitalVector &phi)
        : XCOperator(1, func, phi) { }
    virtual ~XCPotential() { } 

    virtual void setup(double prec);
    virtual void clear();

protected:
    void calcPotential();

    void calcPotentialLDA(int spin);
    void calcPotentialGGA(int spin);
    FunctionTree<3>* calcPotentialGGA(FunctionTreeVector<3> &xc_funcs,
                                      FunctionTreeVector<3> &dRho_a,
                                      FunctionTreeVector<3> &dRho_b);
};

#endif // XCPOTENTIAL_H

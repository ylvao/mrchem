#ifndef DFT_H
#define DFT_H

#include "FockOperator.h"

class DFT : public FockOperator {
public:
    DFT(const MultiResolutionAnalysis<3> &mra,
        KineticOperator &t,
        NuclearPotential &v,
        CoulombOperator &j,
        XCOperator &xc,
        ExchangeOperator *k = 0)
        : FockOperator(mra, &t, &v, &j, k, &xc) { }
    virtual ~DFT() { }
};

#endif // DFT_H

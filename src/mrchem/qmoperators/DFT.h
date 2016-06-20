#ifndef DFT_H
#define DFT_H

#include "FockOperator.h"

class DFT : public FockOperator {
public:
    DFT(KineticOperator &t, NuclearPotential &v,
        CoulombOperator &j, XCOperator &xc, ExchangeOperator *k = 0)
            : FockOperator(&t, &v, &j, k, &xc) { }
    virtual ~DFT() { }
};

#endif // DFT_H

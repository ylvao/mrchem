#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include "FockOperator.h"

class HartreeFock : public FockOperator {
public:
    HartreeFock(KineticOperator &t,
                NuclearPotential &v,
                CoulombOperator &j,
                ExchangeOperator &k)
            : FockOperator(&t, &v, &j, &k) { }
    virtual ~HartreeFock() { }
};

#endif // HARTREEFOCK_H

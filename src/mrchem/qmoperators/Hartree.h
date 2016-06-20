#ifndef HARTREE_H
#define HARTREE_H

#include "FockOperator.h"

class Hartree : public FockOperator {
public:
    Hartree(KineticOperator &t, NuclearPotential &v, CoulombOperator &j) 
        : FockOperator(&t, &v, &j) { }
    virtual ~Hartree() { }
};

#endif // HARTREE_H

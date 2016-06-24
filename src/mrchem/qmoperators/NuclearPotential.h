#ifndef NUCLEARPOTENTIAL_H
#define NUCLEARPOTENTIAL_H

#include "Potential.h"
#include "NuclearFunction.h"
#include "GridCleaner.h"
#include "MWProjector.h"

class NuclearPotential : public Potential {
public:
    NuclearPotential(double build_prec,
                     const MultiResolutionAnalysis<3> &mra,
                     Nuclei &nucs);
    NuclearPotential &operator=(const NuclearPotential &pot) { NOT_IMPLEMENTED_ABORT; }
    virtual ~NuclearPotential();

    virtual void setup(double prec);
    virtual void clear();

    using Potential::operator();
    using Potential::adjoint;

protected:
    NuclearFunction nuc_func;
    GridCleaner<3> clean;
    MWProjector<3> project;
};

#endif // NUCLEARPOTENTIAL_H

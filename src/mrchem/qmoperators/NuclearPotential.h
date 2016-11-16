#ifndef NUCLEARPOTENTIAL_H
#define NUCLEARPOTENTIAL_H

#include "Potential.h"
#include "NuclearFunction.h"
#include "MWProjector.h"

class NuclearPotential : public Potential {
public:
    NuclearPotential(double prec, Nuclei &nucs) : nuc_func(nucs, prec) { }
    NuclearPotential &operator=(const NuclearPotential &pot) { NOT_IMPLEMENTED_ABORT; }
    virtual ~NuclearPotential() { }

    virtual void setup(double prec) {
        Timer timer;
        Potential::setup(prec);
        if (this->real == 0) {
            this->allocReal();
            MWProjector<3> project(this->apply_prec, this->max_scale);
            project(*this->real, this->nuc_func);
        } else {
            NOT_IMPLEMENTED_ABORT;
        }
        timer.stop();
        int n = getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(0, "Nuclear potential", n, t);
    }
    virtual void clear() { Potential::clear(); }

    using Potential::operator();
    using Potential::adjoint;

protected:
    NuclearFunction nuc_func;
};

#endif // NUCLEARPOTENTIAL_H

#ifndef NUCLEARPOTENTIAL_H
#define NUCLEARPOTENTIAL_H

#include "QMTensorOperator.h"
#include "QMPotential.h"
#include "Nucleus.h"
#include "NuclearFunction.h"
#include "MWProjector.h"

class NuclearPotential : public RankZeroTensorOperator {
public:
    NuclearPotential(double Z, const double *R = 0, double S = 1.0e-7) {
        nuc_func.push_back(Z, R, S);
        initializeTensorOperator();
    }
    NuclearPotential(double prec, const Nuclei &nucs)
            : nuclei(nucs), nuc_func(nucs, prec) {
        initializeTensorOperator();
    }
    virtual ~NuclearPotential() { }

    virtual void setup(double prec) {
        this->nuc_pot.setup(prec);
        projectPotential(this->nuc_pot, this->nuc_func);
    }
    virtual void clear() {
        this->nuc_pot.clear();
    }

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

protected:
    Nuclei nuclei;
    NuclearFunction nuc_func;
    QMPotential nuc_pot;

    void initializeTensorOperator() {
        RankZeroTensorOperator &h = *this;
        h = nuc_pot;
    }

    void projectPotential(QMPotential &V, NuclearFunction &f) {
        V.allocReal();

        Timer timer;
        MWProjector<3> project(V.getApplyPrec(), V.getMaxScale());
        project(V.real(), f);
        timer.stop();

        int n = V.getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(0, "Nuclear potential", n, t);
    }
};

#endif // NUCLEARPOTENTIAL_H

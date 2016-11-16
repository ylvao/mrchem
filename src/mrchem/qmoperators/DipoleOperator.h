#ifndef DIPOLEOPERATOR_H
#define DIPOLEOPERATOR_H

#include "Potential.h"
#include "Nucleus.h"
#include "MWProjector.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class DipoleOperator : public Potential {
public:
    DipoleOperator(int dir, double r_0) {
        if (dir < 0 or dir > 2) MSG_ERROR("Invalid direction");

        this->func = [dir, r_0] (const double *r) -> double {
            return r[dir] - r_0;
        };
    }

    virtual ~DipoleOperator() { }

    virtual void setup(double prec) {
        Potential::setup(prec);

        this->real = new FunctionTree<3>(*MRA);
        this->imag = 0;

        MWProjector<3> project(this->apply_prec, this->max_scale);
        project(*this->real, this->func);
    }

    virtual void clear() { Potential::clear(); }

    double operator() (const Nucleus &nuc) {
        const double *R = nuc.getCoord();
        return this->func(R);
    }

    using Potential::operator();

protected:
    std::function<double (const double *r)> func;
};

#endif // DIPOLEOPERATOR_H


#ifndef DIPOLEOPERATOR_H
#define DIPOLEOPERATOR_H

#include "Potential.h"
#include "Nucleus.h"
#include "MWProjector.h"
#include "mrchem.h"

class DipoleOperator : public Potential {
public:
    DipoleOperator(int dir, double r_0)
            : Potential(),
              project(*MRA, -1.0) {
        if (dir < 0 or dir > 2) MSG_ERROR("Invalid direction");

        this->func = [dir, r_0] (const double *r) -> double {
            return r[dir] - r_0;
        };
    }

    virtual ~DipoleOperator() {
    }

    virtual void setup(double prec) {
        Potential::setup(prec);
        this->project.setPrecision(prec);
        this->real = this->project(this->func);
        this->imag = 0;
    }

    virtual void clear() {
        this->project.setPrecision(-1.0);
        Potential::clear();
    }

    double operator() (const Nucleus &nuc) {
        const double *R = nuc.getCoord();
        return this->func(R);
    }

    using Potential::operator();

protected:
    MWProjector<3> project;
    std::function<double (const double *r)> func;
};

#endif // DIPOLEOPERATOR_H


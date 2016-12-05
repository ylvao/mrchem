#ifndef DIPOLEOPERATOR_H
#define DIPOLEOPERATOR_H

#include "Potential.h"
#include "Nucleus.h"
#include "MWProjector.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class DipoleOperator : public Potential {
public:
    DipoleOperator(int i, const double *o = 0) {
        setCoord(this->r_O, o);
        setFunction(i);
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

    double trace(const Nuclei &nucs) {
        double result = 0.0;
        for (int k = 0; k < nucs.size(); k++) {
            result += trace(nucs[k]);
        }
        return result;
    }
    double trace(const Nucleus &nuc) {
        double Z = nuc.getCharge();
        const double *R = nuc.getCoord();
        return -Z*this->func(R);
    }

    using Potential::operator();
    using Potential::adjoint;
    using Potential::trace;

protected:
    double r_O[3];
    std::function<double (const double *r)> func;

    void setCoord(double *out, const double *inp) {
        if (inp != 0) {
            out[0] = inp[0];
            out[1] = inp[1];
            out[2] = inp[2];
        } else {
            out[0] = 0.0;
            out[1] = 0.0;
            out[2] = 0.0;
        }
    }

    void setFunction(int i) {
        const double c = -1.0;
        const double *o = this->r_O;
        switch (i) {
        case 0:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return x;
            };
            break;
        case 1:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return y;
            };
            break;
        case 2:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return z;
            };
            break;
        default:
            MSG_FATAL("Invalid Cartesian component");
        };
    }
};

typedef class DipoleOperator PositionOperator;

#endif // DIPOLEOPERATOR_H


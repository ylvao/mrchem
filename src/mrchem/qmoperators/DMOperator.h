#ifndef DMOPERATOR_H
#define DMOPERATOR_H

#include "Potential.h"
#include "MWProjector.h"

class DMOperator : public Potential {
public:
    DMOperator(int i, int j, const double *o = 0) {
        setOrigin(o);
        setFunction(i, j);
    }
    virtual ~DMOperator() { }

    virtual void setup(double prec) {
        Potential::setup(prec);

        this->real = new FunctionTree<3>(*MRA);
        this->imag = 0;

        MWProjector<3> project(this->apply_prec, this->max_scale);
        project(*this->real, this->func);
    }
    virtual void clear() { Potential::clear(); }

    using Potential::operator();

protected:
    double r_O[3];
    std::function<double (const double *r)> func;

    void setOrigin(const double *o) {
        if (o != 0) {
            this->r_O[0] = o[0];
            this->r_O[1] = o[1];
            this->r_O[2] = o[2];
        } else {
            this->r_O[0] = 0;
            this->r_O[1] = 0;
            this->r_O[2] = 0;
        }
    }

    void setFunction(int i, int j) {
        const double c = -1.0/4.0;
        const double *o = this->r_O;
        int ij = 3*i + j;
        switch (ij) {
        case 0:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return c*(y*y + z*z);
            };
            break;
        case 1:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(x*y);
            };
            break;
        case 2:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(x*z);
            };
            break;
        case 3:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(y*x);
            };
            break;
        case 4:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return c*(x*x + z*z);
            };
            break;
        case 5:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(y*z);
            };
            break;
        case 6:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(z*x);
            };
            break;
        case 7:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return -c*(z*y);
            };
            break;
        case 8:
            this->func = [c, o] (const double *r) -> double {
                double x=r[0]-o[0], y=r[1]-o[1], z=r[2]-o[2];
                return c*(x*x + y*y);
            };
            break;
        default:
            MSG_FATAL("Invalid Cartesian component");
        };
    }
};

#endif // DMOPERATOR_H


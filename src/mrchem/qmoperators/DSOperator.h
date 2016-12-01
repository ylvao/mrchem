#ifndef DSOPERATOR_H
#define DSOPERATOR_H

#include "Potential.h"
#include "NuclearFunction.h"

class DSOperator : public Potential {
public:
    DSOperator(int i, int j, const double *o = 0, const double *k = 0) {
        setCoord(this->r_O, o);
        setCoord(this->r_K, k);
        setFunction(i, j);
    }
    virtual ~DSOperator() { }

    void setup(double prec) {
        Potential::setup(prec);

        this->real = new FunctionTree<3>(*MRA);
        this->imag = 0;

        MWProjector<3> project(this->apply_prec, this->max_scale);
        project(*this->real, this->func);
    }
    virtual void clear() { Potential::clear(); }

    using Potential::operator();
    using Potential::adjoint;

protected:
    double r_O[3];
    double r_K[3];
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

    void setFunction(int i, int j) {
        NuclearFunction r_m1; // 1.0/|r - r_K|
        r_m1.push_back(1.0, this->r_K, 1.0e-7);

        const double mAlpha = 7.2973525664;
        const double c = -pow(mAlpha, 2.0)/2.0;
        const double *o = this->r_O;
        const double *k = this->r_K;
        int ij = 3*i + j;
        switch (ij) {
        case 0:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return c*(y_o*y_k + z_o*z_k)*r_m3;
            };
            break;
        case 1:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(x_o*y_k)*r_m3;
            };
            break;
        case 2:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(x_o*z_k)*r_m3;
            };
            break;
        case 3:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(y_o*x_k)*r_m3;
            };
            break;
        case 4:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return c*(x_o*x_k + z_o*z_k)*r_m3;
            };
            break;
        case 5:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(y_o*z_k)*r_m3;
            };
            break;
        case 6:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(z_o*x_k)*r_m3;
            };
            break;
        case 7:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return -c*(z_o*y_k)*r_m3;
            };
            break;
        case 8:
            this->func = [c, o, k, r_m1] (const double *r) -> double {
                double x_o=r[0]-o[0], y_o=r[1]-o[1], z_o=r[2]-o[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double r_m3 = pow(r_m1.evalf(r), 3.0); 
                return c*(x_o*x_k + y_o*y_k)*r_m3;
            };
            break;
        default:
            MSG_FATAL("Invalid Cartesian component");
        };
    }
};

#endif // DSOPERATOR_H

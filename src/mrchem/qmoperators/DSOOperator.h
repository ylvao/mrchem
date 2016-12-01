#ifndef DSOOPERATOR_H
#define DSOOPERATOR_H

#include "Potential.h"
#include "NuclearFunction.h"

class DSOOperator : public Potential {
public:
    DSOOperator(int i, int j, const double *k = 0, const double *l = 0 ) {
        setCoord(this->r_K, k);
        setCoord(this->r_L, l);
        setFunction(i, j);
    }
    virtual ~DSOOperator() { }

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
    double r_K[3];
    double r_L[3];
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
        NuclearFunction K_m1; // 1.0/|r - r_K|
        NuclearFunction L_m1; // 1.0/|r - r_L|
        K_m1.push_back(1.0, this->r_K, 1.0e-7);
        L_m1.push_back(1.0, this->r_L, 1.0e-7);

        const double mAlpha = 7.2973525664;
        const double c = -pow(mAlpha, 4.0)/2.0;
        const double *k = this->r_K;
        const double *l = this->r_L;
        int ij = 3*i + j;
        switch (ij) {
        case 0:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return c*(y_k*y_l + z_k*z_l)*K_m3*L_m3;
            };
            break;
        case 1:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(x_k*y_l)*K_m3*L_m3;
            };
            break;
        case 2:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(x_k*z_l)*K_m3*L_m3;
            };
            break;
        case 3:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(y_k*x_l)*K_m3*L_m3;
            };
            break;
        case 4:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return c*(x_k*x_l + z_k*z_l)*K_m3*L_m3;
            };
            break;
        case 5:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(y_k*z_l)*K_m3*L_m3;
            };
            break;
        case 6:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(z_k*x_l)*K_m3*L_m3;
            };
            break;
        case 7:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return -c*(z_k*y_l)*K_m3*L_m3;
            };
            break;
        case 8:
            this->func = [c, k, l, K_m1, L_m1] (const double *r) -> double {
                double x_l=r[0]-l[0], y_l=r[1]-l[1], z_l=r[2]-l[2];
                double x_k=r[0]-k[0], y_k=r[1]-k[1], z_k=r[2]-k[2];
                double K_m3 = pow(K_m1.evalf(r), 3.0);
                double L_m3 = pow(L_m1.evalf(r), 3.0);
                return c*(x_k*x_l + y_k*y_l)*K_m3*L_m3;
            };
            break;
        default:
            MSG_FATAL("Invalid Cartesian component");
        };
    }
};

#endif // DIASHIELDOPERATOR_H

#ifndef H_M_PSO_H
#define H_M_PSO_H

#include "QMTensorOperator.h"
#include "AngularMomentumOperator.h"
#include "NuclearPotential.h"

class H_M_pso : public RankOneTensorOperator<3> {
public:
    H_M_pso(DerivativeOperator<3> &d, const double *k = 0) : l(d, k), r_m1(1.0, k) {
        initializeTensorOperator();
    }
    virtual ~H_M_pso() { }

    void setup(double prec) {
        this->l.setup(prec);
        this->r_m1.setup(prec);
    }
    void clear() {
        this->l.clear();
        this->r_m1.clear();
    }

protected:
    NuclearPotential r_m1;
    AngularMomentumOperator l;

    void initializeTensorOperator() {
        static double alpha = 7.2973525664;
        RankZeroTensorOperator &l_x = this->l[0];
        RankZeroTensorOperator &l_y = this->l[1];
        RankZeroTensorOperator &l_z = this->l[2];
        RankZeroTensorOperator r_m3 = r_m1*r_m1*r_m1;

        RankOneTensorOperator<3> &h = (*this);
        h[0] = (alpha*alpha)*r_m3*l_x;
        h[1] = (alpha*alpha)*r_m3*l_y;
        h[2] = (alpha*alpha)*r_m3*l_z;
    }
};

#endif // H_M_PSO_H


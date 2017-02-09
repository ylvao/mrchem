#ifndef H_E_DIP_H
#define H_E_DIP_H

#include "QMTensorOperator.h"
#include "PositionOperator.h"

class H_E_dip : public RankOneTensorOperator<3> {
public:
    H_E_dip(const double *o = 0) : r_x(0, o), r_y(1, o), r_z(2, o) {
        initializeTensorOperator();
    }
    virtual ~H_E_dip() { }

    virtual void setup(double prec) {
        this->r_x.setup(prec);
        this->r_y.setup(prec);
        this->r_z.setup(prec);
    }

    virtual void clear() {
        this->r_x.clear();
        this->r_y.clear();
        this->r_z.clear();
    }

    Eigen::VectorXd trace(const Nuclei &nucs) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
        for (int k = 0; k < nucs.size(); k++) {
            result += trace(nucs[k]);
        }
        return result;
    }
    Eigen::VectorXd trace(const Nucleus &nuc) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(3);
        double Z = nuc.getCharge();
        const double *R = nuc.getCoord();
        result(0) = -Z*this->r_x.func(R);
        result(1) = -Z*this->r_y.func(R);
        result(2) = -Z*this->r_z.func(R);
        return result;
    }

    using RankOneTensorOperator<3>::trace;

protected:
    PositionOperator r_x;
    PositionOperator r_y;
    PositionOperator r_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_x;
        h[1] = r_y;
        h[2] = r_z;
    }
};

#endif // H_E_DIP_H


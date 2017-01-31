#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMTensorOperator.h"
#include "QMDerivative.h"

class MomentumOperator : public RankOneTensorOperator<3> {
public:
    MomentumOperator() : d_x(0), d_y(1), d_z(2) {
        initializeTensorOperator();
    }
    virtual ~MomentumOperator() { }

    virtual void setup(double prec) {
        this->d_x.setup(prec);
        this->d_y.setup(prec);
        this->d_z.setup(prec);
    }
    virtual void clear() {
        this->d_x.clear();
        this->d_y.clear();
        this->d_z.clear();
    }

protected:
    QMDerivative d_x;
    QMDerivative d_y;
    QMDerivative d_z;

    void initializeTensorOperator() {
        std::complex<double> i(0.0, 1.0);

        RankOneTensorOperator<3> &h = *this;
        h[0] = i*d_x;
        h[1] = i*d_y;
        h[2] = i*d_z;
    }
};

#endif // MOMENTUMOPERATOR_H


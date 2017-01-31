#ifndef POSITIONOPERATOR_H
#define POSITIONOPERATOR_H

#include "QMTensorOperator.h"
#include "QMPotential.h"
#include "MWProjector.h"

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const double *o = 0) {
        setCoord(this->r_O, o);
        initializeFunctions();
        initializeTensorOperator();
    }
    virtual ~PositionOperator() { }

    virtual void setup(double prec) {
        this->r_x.setup(prec);
        this->r_y.setup(prec);
        this->r_z.setup(prec);

        projectPotential(this->r_x, this->f_x);
        projectPotential(this->r_y, this->f_y);
        projectPotential(this->r_z, this->f_z);
    }

    virtual void clear() {
        this->r_x.clear();
        this->r_y.clear();
        this->r_z.clear();
    }

protected:
    double r_O[3];
    QMPotential r_x;
    QMPotential r_y;
    QMPotential r_z;
    std::function<double (const double *r)> f_x;
    std::function<double (const double *r)> f_y;
    std::function<double (const double *r)> f_z;

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

    void initializeFunctions() {
        const double *o = this->r_O;
        this->f_x = [o] (const double *r) -> double { return r[0]-o[0]; };
        this->f_y = [o] (const double *r) -> double { return r[1]-o[1]; };
        this->f_z = [o] (const double *r) -> double { return r[2]-o[2]; };
    }

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = *this;
        h[0] = r_x;
        h[1] = r_y;
        h[2] = r_z;
    }

    void projectPotential(QMPotential &V, std::function<double (const double *r)> &f) {
        V.allocReal();

        Timer timer;
        MWProjector<3> project(V.getApplyPrec(), V.getMaxScale());
        project(V.real(), f);
        timer.stop();

        int n = V.getNNodes();
        double t = timer.getWallTime();
        TelePrompter::printTree(1, "Position operator", n, t);
    }
};

#endif // POSITIONOPERATOR_H


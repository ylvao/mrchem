#ifndef POSITIONOPERATOR_H
#define POSITIONOPERATOR_H

#include "QMPotential.h"
#include "QMTensorOperator.h"


class QMPosition : public QMPotential {
public:
    QMPosition(int d, const double *o = 0) {
        setPosition(this->r_O, o);
        setFunction(d, this->r_O);
    }
    virtual ~QMPosition() { }

    virtual void setup(double prec);
    virtual void clear();

protected:
    double r_O[3];
    std::function<double (const double *r)> func;

    void setPosition(double *out, const double *inp);
    void setFunction(int d, const double *o);
};

class QMPositionX : public QMPosition {
public:
    QMPositionX(const double *o = 0) : QMPosition(0, o) { }
    virtual ~QMPositionX() { }
};

class QMPositionY : public QMPosition {
public:
    QMPositionY(const double *o = 0) : QMPosition(1, o) { }
    virtual ~QMPositionY() { }
};

class QMPositionZ : public QMPosition {
public:
    QMPositionZ(const double *o = 0) : QMPosition(2, o) { }
    virtual ~QMPositionZ() { }
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const double *o = 0) : r_x(o), r_y(o), r_z(o) {
        initializeTensorOperator();
    }
    virtual ~PositionOperator() { }

protected:
    QMPositionX r_x;
    QMPositionY r_y;
    QMPositionZ r_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = r_x;
        h[1] = r_y;
        h[2] = r_z;
    }
};

#endif // POSITIONOPERATOR_H

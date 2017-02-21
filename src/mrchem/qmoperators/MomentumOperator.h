#ifndef MOMENTUMOPERATOR_H
#define MOMENTUMOPERATOR_H

#include "QMDerivative.h"
#include "QMTensorOperator.h"

class QMMomentum : public QMDerivative {
public:
    QMMomentum(int d, DerivativeOperator<3> &D) : QMDerivative(d, D) { }
    virtual ~QMMomentum() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

protected:
    virtual bool isReal() const { return false; }
    virtual bool isImag() const { return true; }
};

class QMMomentumX : public QMMomentum {
public:
    QMMomentumX(DerivativeOperator<3> &D) : QMMomentum(0, D) { }
    virtual ~QMMomentumX() { }
};

class QMMomentumY : public QMMomentum {
public:
    QMMomentumY(DerivativeOperator<3> &D) : QMMomentum(1, D) { }
    virtual ~QMMomentumY() { }
};

class QMMomentumZ : public QMMomentum {
public:
    QMMomentumZ(DerivativeOperator<3> &D) : QMMomentum(2, D) { }
    virtual ~QMMomentumZ() { }
};

class MomentumOperator : public RankOneTensorOperator<3> {
public:
    MomentumOperator(DerivativeOperator<3> &D) : p_x(D), p_y(D), p_z(D) {
        initializeTensorOperator();
    }
    virtual ~MomentumOperator() { }

protected:
    QMMomentumX p_x;
    QMMomentumY p_y;
    QMMomentumZ p_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = p_x;
        h[1] = p_y;
        h[2] = p_z;
    }
};

#endif // MOMENTUMOPERATOR_H


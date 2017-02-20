#ifndef SPINOPERATOR_H
#define SPINOPERATOR_H

#include "QMOperator.h"
#include "QMTensorOperator.h"

extern MultiResolutionAnalysis<3> *MRA; // Global MRA

class QMSpinX : public QMOperator {
public:
    QMSpinX() : QMOperator(MRA->getMaxScale()) { }
    virtual ~QMSpinX() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class QMSpinY : public QMOperator {
public:
    QMSpinY() : QMOperator(MRA->getMaxScale()) { }
    virtual ~QMSpinY() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class QMSpinZ : public QMOperator {
public:
    QMSpinZ() : QMOperator(MRA->getMaxScale()) { }
    virtual ~QMSpinZ() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital* operator()(Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    using QMOperator::operator();
    using QMOperator::adjoint;
};

class SpinOperator : public RankOneTensorOperator<3> {
public:
    SpinOperator() { initializeTensorOperator(); }
    virtual ~SpinOperator() { }

protected:
    QMSpinX s_x;
    QMSpinY s_y;
    QMSpinZ s_z;

    void initializeTensorOperator() {
        RankOneTensorOperator<3> &h = (*this);
        h[0] = s_x;
        h[1] = s_y;
        h[2] = s_z;
    }
};

#endif // SPINOPERATOR_H

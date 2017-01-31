#ifndef QMOPERATOR_H
#define QMOPERATOR_H

#include <Eigen/Core>

#include "TelePrompter.h"
#include "constants.h"

class Orbital;

/** \brief Quantum mechanical Hermitian operators

    Base class to handle operators and their application in the QM sense

  */
class QMOperator {
public:
    QMOperator(int ms) : max_scale(ms), apply_prec(-1.0) { }
    QMOperator(const QMOperator &oper) : max_scale(oper.max_scale), apply_prec(oper.apply_prec) { }
    QMOperator& operator=(const QMOperator &inp) { this->apply_prec = inp.apply_prec; return *this; }
    virtual ~QMOperator() { }

    int getMaxScale() const { return this->max_scale; }
    double getApplyPrec() const { return this->apply_prec; }

    virtual void setup(double prec) {
        if (this->apply_prec < 0.0) { 
            this->apply_prec = prec;
        } else if (fabs(prec - this->apply_prec) > MachineZero) {
            MSG_ERROR("Clear operator before setup with different prec!");
        }
    }
    virtual void clear() { this->apply_prec = -1.0; }

    virtual Orbital* operator() (Orbital &phi) = 0;
    virtual Orbital* adjoint(Orbital &phi) = 0;

protected:
    const int max_scale;
    double apply_prec;
};

#endif // QMOPERATOR_H


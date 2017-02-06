#ifndef QMOPERATOR_H
#define QMOPERATOR_H

#include <Eigen/Core>

#include "TelePrompter.h"
#include "constants.h"

class Orbital;
class OrbitalVector;

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

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &phi) = 0;
    virtual Orbital* adjoint(Orbital &phi) = 0;

    virtual double operator() (Orbital &phi_i, Orbital &phi_j);
    virtual double adjoint(Orbital &phi_i, Orbital &phi_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    const int max_scale;
    double apply_prec;
};

#endif // QMOPERATOR_H


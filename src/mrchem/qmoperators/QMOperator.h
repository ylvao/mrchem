#ifndef QMOPERATOR_H
#define QMOPERATOR_H

#include <Eigen/Core>

class OrbitalVector;
class Orbital;

/** \brief Quantum mechanical Hermitian operators

    Base class to handle operators and their application in the QM sense

  */
class QMOperator {
public:
    QMOperator(int ms) : max_scale(ms), apply_prec(-1.0) { }
    QMOperator(const QMOperator &oper) : max_scale(oper.max_scale), apply_prec(oper.apply_prec) { }
    virtual ~QMOperator() { }

    virtual void setup(double prec) { this->apply_prec = prec; }
    virtual void clear() { this->apply_prec = -1.0; }

    virtual Orbital* operator() (Orbital &phi_p) = 0;
    virtual Orbital* adjoint(Orbital &phi_p) = 0;

    virtual double operator() (Orbital &phi_i, Orbital &phi_j);
    virtual double adjoint(Orbital &phi_i, Orbital &phi_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    double trace(Orbital &phi_p, Orbital *x_p, Orbital *y_p);
    double trace(OrbitalVector &phi, OrbitalVector *x, OrbitalVector *y);
    double trace(OrbitalVector &phi);

protected:
    const int max_scale;
    double apply_prec;
};

#endif // QMOPERATOR_H


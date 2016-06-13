#ifndef QMOPERATOR_H
#define QMOPERATOR_H

#include <Eigen/Core>
#include <string>

#include "TelePrompter.h"

class OrbitalVector;
class Orbital;

/** \brief Quantum mechanical operators

    Base class to handle operators and their application in the QM sense

  */
class QMOperator {
public:
    QMOperator(bool im = false)
            : imaginary(im),
              precision(-1.0) {
    }
    virtual ~QMOperator() { }

    virtual void rotate(double prec, Eigen::MatrixXd &U) { }

    virtual void setupUnperturbed(double prec = -1.0) { }
    virtual void clearUnperturbed() { }

    virtual void setup(double prec = -1.0) { this->precision = prec; }
    virtual void clear() { this->precision = -1.0; }

    bool isImaginary() const { return this->imaginary; }

    virtual int printTreeSizes() const;

    virtual Orbital* operator() (Orbital &orb_p) = 0;
    virtual Orbital* adjoint(Orbital &orb_p) = 0;

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    double calcProperty(Orbital &phi_p, Orbital *x_p, Orbital *y_p);
    double calcProperty(OrbitalVector &phi, OrbitalVector *x, OrbitalVector *y);

protected:
    bool imaginary;
    double precision;
};

#endif // QMOPERATOR_H


#ifndef QMOPERATOR_H
#define QMOPERATOR_H

#include <Eigen/Core>
#include <string>

#include "GridCleaner.h"
#include "GridGenerator.h"

class OrbitalVector;
class Orbital;

/** \brief Quantum mechanical Hermitian operators

    Base class to handle operators and their application in the QM sense

  */
class QMOperator {
public:
    QMOperator() : apply_prec(-1.0), clean(-1.0) { }
    virtual ~QMOperator() { }

    virtual void rotate(Eigen::MatrixXd &U) { }

    virtual void setup(double prec);
    virtual void clear();

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
    double apply_prec;
    GridCleaner<3> clean;
    GridGenerator<3> grid;
};

#endif // QMOPERATOR_H


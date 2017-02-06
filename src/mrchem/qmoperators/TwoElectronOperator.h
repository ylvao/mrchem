#ifndef TWOELECTRONOPERATOR_H
#define TWOELECTRONOPERATOR_H

#include <Eigen/Core>

#include "QMOperator.h"

class OrbitalVector;

class TwoElectronOperator : public QMOperator {
public:
    TwoElectronOperator(int ms, OrbitalVector &phi)
        : QMOperator(ms), orbitals(&phi) { }
    virtual ~TwoElectronOperator() { }

    virtual void rotate(Eigen::MatrixXd &U) { NOT_REACHED_ABORT; }

protected:
    OrbitalVector *orbitals;
};

#endif // TWOELECTRONOPERATOR_H



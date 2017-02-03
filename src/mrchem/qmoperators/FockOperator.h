#ifndef FOCKOPERATOR_H
#define FOCKOPERATOR_H

#include <vector>

#include "QMTensorOperator.h"

class QMOperatorExp;
class NuclearPotential;
class KineticOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class Orbital;
class SCFEnergy;

class FockOperator : public RankZeroTensorOperator {
public:
    FockOperator(KineticOperator *t = 0,
                 NuclearPotential *v = 0,
                 CoulombOperator *j = 0,
                 ExchangeOperator *k = 0,
                 XCOperator *xc = 0);
    virtual ~FockOperator();

    void addPerturbationOperator(QMOperatorExp &h_1) { this->H_1 = &h_1; }

    KineticOperator *getKineticOperator() { return this->T; }
    NuclearPotential *getNuclearPotential() { return this->V; }
    CoulombOperator *getCoulombOperator() { return this->J; }
    ExchangeOperator *getExchangeOperator() { return this->K; }
    XCOperator *getXCOperator() { return this->XC; }
    QMOperatorExp *getPerturbationOperator() { return this->H_1; }

    virtual void rotate(Eigen::MatrixXd &U);

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    SCFEnergy trace(OrbitalVector &phi, Eigen::MatrixXd &F);

    Orbital* applyKinetic(Orbital &orb_p);
    double applyKinetic(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyKinetic(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    Orbital* applyAdjointKinetic(Orbital &orb_p);
    double applyAdjointKinetic(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyAdjointKinetic(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    Orbital* applyPotential(Orbital &orb_p);
    double applyPotential(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    Orbital* applyAdjointPotential(Orbital &orb_p);
    double applyAdjointPotential(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyAdjointPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    KineticOperator *T;
    NuclearPotential *V;
    CoulombOperator *J;
    ExchangeOperator *K;
    XCOperator *XC;
    QMOperatorExp *H_1;   // First order perturbation operators
};

#endif // FOCKOPERATOR_H


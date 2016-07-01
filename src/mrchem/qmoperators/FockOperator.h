#ifndef FOCKOPERATOR_H
#define FOCKOPERATOR_H

#include <vector>

#include "QMOperator.h"
#include "OrbitalAdder.h"

class Molecule;
class Orbital;
class NuclearPotential;
class KineticOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;

class FockOperator : public QMOperator {
public:
    FockOperator(const MultiResolutionAnalysis<3> &mra,
                 KineticOperator *t = 0,
                 NuclearPotential *v = 0,
                 CoulombOperator *j = 0,
                 ExchangeOperator *k = 0,
                 XCOperator *xc = 0);
    virtual ~FockOperator();

    int getNPerturbations() const { return this->H_1.size(); }
    void addPerturbationOperator(QMOperator &h_1) { this->H_1.push_back(&h_1); }
    QMOperator& getPerturbationOperator(int i) { return *this->H_1[i]; }
    const QMOperator& getPerturbationOperator(int i) const { return *this->H_1[i]; }

    KineticOperator *getKineticOperator() { return this->T; }
    NuclearPotential *getNuclearPotential() { return this->V; }
    CoulombOperator *getCoulombOperator() { return this->J; }
    ExchangeOperator *getExchangeOperator() { return this->K; }
    XCOperator *getXCOperator() { return this->XC; }

    int printTreeSizes() const;

    virtual void rotate(Eigen::MatrixXd &U);

    virtual void setup(double prec);
    virtual void clear();

    virtual Orbital* operator() (Orbital &orb_p);
    virtual Orbital* adjoint(Orbital &orb_p);

    virtual double operator() (Orbital &orb_i, Orbital &orb_j);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

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

    Orbital* applyPerturbations(Orbital &orb_p);
    double applyPerturbations(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyPerturbations(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

    Orbital* applyAdjointPerturbations(Orbital &orb_p);
    double applyAdjointPerturbations(Orbital &orb_i, Orbital &orb_j);
    Eigen::MatrixXd applyAdjointPerturbations(OrbitalVector &i_orbs, OrbitalVector &j_orbs);

protected:
    OrbitalAdder add;
    KineticOperator *T;
    NuclearPotential *V;
    CoulombOperator *J;
    ExchangeOperator *K;
    XCOperator *XC;
    std::vector<QMOperator *> H_1;   // First order perturbation operators
};

#endif // FOCKOPERATOR_H


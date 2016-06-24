#ifndef GROUNDSTATESOLVER_H
#define GROUNDSTATESOLVER_H

#include <Eigen/Core>

#include "SCF.h"
#include "SCFEnergy.h"

class KAIN;

class GroundStateSolver : public SCF {
public:
    GroundStateSolver(const MultiResolutionAnalysis<3> &mra,
                      HelmholtzOperatorSet &h);
    virtual ~GroundStateSolver();

protected:
    std::vector<SCFEnergy> energy;

    FockOperator *fOper_n;
    Eigen::MatrixXd *fMat_n;

    OrbitalVector *orbitals_n;
    OrbitalVector *orbitals_np1;
    OrbitalVector *dOrbitals_n;

    virtual Orbital* getHelmholtzArgument(int i,
                                          Eigen::MatrixXd &F,
                                          OrbitalVector &phi,
                                          bool adjoint);

    void printProperty() const;
    double calcProperty();
    double calcPropertyError() const;

    void localize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    void diagonalize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    void orthonormalize(FockOperator &fock, Eigen::MatrixXd &F, OrbitalVector &phi);
    Eigen::MatrixXd calcOrthonormalizationMatrix(OrbitalVector &phi);
};

#endif // GROUNDSTATESOLVER_H


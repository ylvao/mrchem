#ifndef ENERGYOPTIMIZER_H
#define ENERGYOPTIMIZER_H

#include "GroundStateSolver.h"

class EnergyOptimizer : public GroundStateSolver {
public:
    EnergyOptimizer(const MultiResolutionAnalysis<3> &mra,
                    HelmholtzOperatorSet &h);
    virtual ~EnergyOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, Eigen::MatrixXd &F,
               FockOperator &fock_np1, OrbitalVector &phi_np1);
    void clear();

    virtual bool optimize();

protected:
    FockOperator *fOper_np1;

    Eigen::MatrixXd calcFockMatrixUpdate();
};

#endif // ENERGYOPTIMIZER_H

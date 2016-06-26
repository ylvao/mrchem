#ifndef ORBITALOPTIMIZER_H
#define ORBITALOPTIMIZER_H

#include "GroundStateSolver.h"

class Accelerator;

class OrbitalOptimizer : public GroundStateSolver {
public:
    OrbitalOptimizer(const MultiResolutionAnalysis<3> &mra,
                     HelmholtzOperatorSet &h,
                     Accelerator *k = 0);
    virtual ~OrbitalOptimizer();

    void setup(FockOperator &fock, OrbitalVector &phi, Eigen::MatrixXd &F);
    void clear();

    virtual bool optimize();

protected:
    Accelerator *kain; // Pointer to external object, do not delete!

    void printTreeSizes() const;
};

#endif // ORBITALOPTIMIZER_H

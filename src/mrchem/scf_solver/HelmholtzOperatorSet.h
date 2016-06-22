#ifndef HELMHOLTZOPERATORSET_H
#define HELMHOLTZOPERATORSET_H

#include <Eigen/Core>
#include <vector>

#include "HelmholtzOperator.h"

class HelmholtzOperatorSet {
public:
    HelmholtzOperatorSet(double build,
                         const MultiResolutionAnalysis<3> &mra,
                         double thrs = -1.0);
    virtual ~HelmholtzOperatorSet() { clear(); }

    void initialize(const Eigen::VectorXd &energies);
    void clear();

    int printTreeSizes() const;

    void setThreshold(double thrs) { this->threshold = thrs; }
    double getThreshold() const { return this->threshold; }

    double getLambda(int i) const { return this->lambda[i]; }
    Eigen::VectorXd getLambda() const;
    HelmholtzOperator &getOperator(int i);

private:
    double threshold; //For re-using operators. Negative means always recreate
    double build_prec;
    const MultiResolutionAnalysis<3> MRA;

    std::vector<int> operIdx;
    std::vector<double> lambda;
    std::vector<HelmholtzOperator *> operators;

    int initHelmholtzOperator(double energy);
    void clearUnused();
};


#endif // ORBITALSET_H

#ifndef DENSEXP_H
#define DENSEXP_H

#include <Eigen/Core>
#include <vector>

#include "OrbExp.h"
#include "RepresentableFunction.h"

class DensExp : public RepresentableFunction<3>, public OrbExp {
public:
    DensExp(Intgrl &intgrl, Eigen::MatrixXd &D);
    virtual ~DensExp() { }

    Eigen::MatrixXd &getDensityMatrix() { return this->densityMatrix; }
    GaussExp<3> getAODensExpansion();

    double evalf(const double *r) const;
    void calcScreening(double nStdDev);
    void setScreen(bool screen);

    //bool checkSeedNode(MWNode<3> &node);
    //void calcWaveletCoefs(MWNode<3> &node);

protected:
    Eigen::MatrixXd densityMatrix;

    //int sortDensTerms(MWNode<3> &node,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &outDens,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &inpVals,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &outVals);
};

#endif // DENSEXP_H

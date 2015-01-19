#ifndef DENSEXP_H
#define DENSEXP_H

#include <Eigen/Core>

#include "RepresentableFunction.h"

template<int D> class GaussExp;

class DensExp : public RepresentableFunction<3> {
public:
    DensExp(Intgrl &intgrl, Eigen::MatrixXd &D);
    virtual ~DensExp() { }

    Eigen::MatrixXd &getDensityMatrix() { return &this->densMat; }
    GaussExp<3> getAODensExpansion();
    GaussExp<3> &getAO(int n) { return *this->orbitals[n]; }

    int size() const { return this->orbitals.size(); }
    int getAngularMomentum(int n) const;

    double evalf(const double *r) const;
    void calcScreening(double nStdDev);
    void setScreen(bool screen);

    void rotate(Eigen::MatrixXd &U);

    //bool checkSeedNode(MWNode<3> &node);
    //void calcWaveletCoefs(MWNode<3> &node);

protected:
    bool cartesian;
    Eigen::MatrixXd &densMat;
    std::vector<GaussExp<3> *> orbitals;

    void readAOExpansion(Intgrl &intgrl);
    void transformToSpherical();

    //int sortDensTerms(MWNode<3> &node,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &outDens,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &inpVals,
    //                  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
    //                  Eigen::RowMajor> &outVals);
};

#endif // DENSEXP_H

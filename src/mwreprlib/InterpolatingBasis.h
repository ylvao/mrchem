/*
 *
 *
 *  \date June 2, 2010
 *  \author Stig Rune Jensen \n
 *          CTCC, University of Tromsø
 *
 */

#ifndef INTERPOLATINGBASIS_H
#define INTERPOLATINGBASIS_H

#include "ScalingBasis.h"

class InterpolatingBasis : public ScalingBasis {
public:
    InterpolatingBasis(int order) :	 ScalingBasis(order) {
//        initScalingBasis();
//        preEvaluate();
    }
    virtual ~InterpolatingBasis() { }

//    Eigen::VectorXd calcScalingCoefs(int axis, const SeparableFunction<1> &func,
//        int n, int l) const;
//    Eigen::VectorXd calcScalingCoefs(int axis, const SeparableFunction<2> &func,
//        int n, int l) const;
//    Eigen::VectorXd calcScalingCoefs(int axis, const SeparableFunction<3> &func,
//        int n, int l) const;

//    void calcScalingCoefs(const SeparableFunction<1> &func, int n, const int *l,
//        Eigen::MatrixXd &cfs) const;
//    void calcScalingCoefs(const SeparableFunction<2> &func, int n, const int *l,
//        Eigen::MatrixXd &cfs) const;
//    void calcScalingCoefs(const SeparableFunction<3> &func, int n, const int *l,
//        Eigen::MatrixXd &cfs) const;
protected:
//    void initScalingBasis();
//    void preEvaluate();
};

#endif // INTERPOLATINGBASIS_H

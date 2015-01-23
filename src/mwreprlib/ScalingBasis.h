/*
 *
 *
 *  \date Oct 15, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef SCALINGBASIS_H_
#define SCALINGBASIS_H_

#include "TelePrompter.h"

class ScalingBasis {
public:
    ScalingBasis(int order) {
        if (order < 1) MSG_ERROR("Invalid scaling order");
        this->scalingOrder = order;
        this->quadratureOrder = order + 1;
//        double lB = 0.0;
//        double uB = 1.0;
//        this->setAllBounds(&lB, &uB);
    }
    virtual ~ScalingBasis() { }

//    double evalf(int k, double x) const {
//        NOT_IMPLEMENTED_ABORT;
//        assert(k >= 0 and k <= this->scalingOrder);
//        return this->getFunc(k).evalf(x);
//    }
//    void evalf(const double *x, Eigen::MatrixXd &vals) const {
//        int dim = vals.cols();
//        for (int i = 0; i < dim; i++) {
//            for (int k = 0; k < this->scalingOrder + 1; k++) {
//                vals(k, i) = this->getFunc(k).evalf(x[i]);
//            }
//        }
//    }

//    virtual Eigen::VectorXd calcScalingCoefs(int axis,
//            const SeparableFunction<1> &func, int n, int l) const = 0;
//    virtual Eigen::VectorXd calcScalingCoefs(int axis,
//            const SeparableFunction<2> &func, int n, int l) const = 0;
//    virtual Eigen::VectorXd calcScalingCoefs(int axis,
//            const SeparableFunction<3> &func, int n, int l) const = 0;

//    virtual void calcScalingCoefs(const SeparableFunction<1> &func, int n,
//            const int *l, Eigen::MatrixXd &cfs) const = 0;
//    virtual void calcScalingCoefs(const SeparableFunction<2> &func, int n,
//            const int *l, Eigen::MatrixXd &cfs) const = 0;
//    virtual void calcScalingCoefs(const SeparableFunction<3> &func, int n,
//            const int *l, Eigen::MatrixXd &cfs) const = 0;

//    double getVal(int k, int i) const { return this->preVals(k, i); }
//    int getScalingOrder() const { return this->scalingOrder; }
//    int getType() const { return this->type; }
//    int getQuadratureOrder() const { return this->quadratureOrder; }
//    void setQuadratureOrder(int order) {
//        if (order < 1) MSG_ERROR("Quadrature order < 1: " << order);
//        this->quadratureOrder = order;
//        preEvaluate();
//    }
protected:
    int scalingOrder;
    int quadratureOrder;

//    Eigen::MatrixXd preVals;

//    virtual void initScalingBasis() = 0;
//    virtual void preEvaluate() = 0;
};

#endif /* SCALINGBASIS_H_ */

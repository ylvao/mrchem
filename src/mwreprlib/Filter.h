/*
 *
 *
 *  \date Jul 31, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif Abstract base class for filters and filter oprtations.
 */

#ifndef FILTER_H_
#define FILTER_H_

#include <Eigen/Core>

#include "MultiException.h"
#include "constants.h"

class Filter {
public:
    Filter(int k, int t) : type(t), order(k) {
        if (this->order < 1) {
            THROW_ERROR("Filter order < 1. Requested order: " << this->order)
        }
        switch (this->type) {
        case (Interpol):
        case (Legendre):
            break;
        default:
            THROW_ERROR("Unknown filter type: " << type)
        }
    }
    virtual ~Filter() { }
    virtual void apply(Eigen::VectorXd &data) const = 0;
    virtual void applyInverse(Eigen::VectorXd &data) const = 0;
    virtual void apply(Eigen::MatrixXd &data) const = 0;
    virtual void applyInverse(Eigen::MatrixXd &data) const = 0;

    int getOrder() const { return this->order; }
    int getType() const { return this->type; }

    virtual const Eigen::MatrixXd &getFilter() const = 0;
    virtual const Eigen::MatrixXd &getSubFilter(int i, int oper = 0) const = 0;
    virtual const Eigen::MatrixXd &getCompressionSubFilter(int i) const = 0;
    virtual const Eigen::MatrixXd &getReconstructionSubFilter(int i) const = 0;
protected:
    int type;
    int order;
    int dim;
};

#endif /* FILTER_H_ */

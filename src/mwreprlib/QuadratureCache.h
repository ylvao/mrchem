/*
 *
 *
 *  \date Jul 26, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#ifndef QUADRATURECACHE_H_
#define QUADRATURECACHE_H_

#include <Eigen/Core>

#include "GaussQuadrature.h"
#include "ObjectCache.h"

#define getQuadratureCache(X)\
    QuadratureCache &X = QuadratureCache::getInstance()

class QuadratureCache: public ObjectCache<GaussQuadrature> {
public:
    static QuadratureCache &getInstance() {
    NOT_IMPLEMENTED_ABORT;
        static QuadratureCache theQuadratureCache;
        return theQuadratureCache;
    }
    virtual void load(int order);
    virtual GaussQuadrature &get(int order);
    virtual const Eigen::VectorXd &getRoots(int i) {
    NOT_IMPLEMENTED_ABORT;
        return get(i).getRoots();
    }
    virtual const Eigen::VectorXd &getWeights(int i) {
    NOT_IMPLEMENTED_ABORT;
        return get(i).getWeights();
    }
    void setIntervals(int i);
    void setBounds(double a, double b);
    int getIntervals() const {
        return intervals;
    }
    double getUpperBound() const {
        return B;
    }
    double getLowerBound() const {
        return A;
    }
private:
    double A;
    double B;
    int intervals;
    QuadratureCache();
    virtual ~QuadratureCache();
    QuadratureCache(QuadratureCache const &qc):
        ObjectCache<GaussQuadrature>(qc) {
    NOT_IMPLEMENTED_ABORT;
    }
    QuadratureCache &operator=(QuadratureCache const&) {
    NOT_IMPLEMENTED_ABORT;
        return *this;
    }
};

#endif /* QUADRATURECACHE_H_ */

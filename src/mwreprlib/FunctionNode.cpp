#include <Eigen/Core>
#include <set>

#include "MWNode.h"
#include "FunctionNode.h"
#include "FunctionTree.h"
#include "QuadratureCache.h"
#include "MathUtils.h"
#include "config.h"
#include "eigen_disable_warnings.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

/** FunctionNode default constructor.
  * Creates an empty node. */
template<int D>
FunctionNode<D>::FunctionNode() : MWNode<D> () {
}

/** Root node constructor. By default root nodes are initialized to
  * represent functions which are constant zero. */
template<int D>
FunctionNode<D>::FunctionNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx)
        : MWNode<D> (t, nIdx) {
}

/** FunctionNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
FunctionNode<D>::FunctionNode(FunctionNode<D> *p, int cIdx)
        : MWNode<D> (p, cIdx) {
}

/** FunctionNode copy constructor.
  * Make a copy of a node and assign it to another parent.
  * Copying coefficients is optional. */
template<int D>
FunctionNode<D>::FunctionNode(const FunctionNode<D> &nd, FunctionNode<D> *p)
        : MWNode<D> (nd, p) {
}

/** FunctionNode copy constructor.
  *
  * Make a detached copy of a node that is not accessible through any tree. The
  * tree is still accessible from the node, as much of the node parameters are
  * in fact stored in the tree. Copying coefficients is optional. Children
  * nodes are NOT copied recursively. */
template<int D>
FunctionNode<D>::FunctionNode(const FunctionNode<D> &nd, FunctionTree<D> *t)
        : MWNode<D> (nd, t) {
}

/** FunctionNode equals operator.
  * Copying the content of a node, not its location. */
template<int D>
FunctionNode<D> &FunctionNode<D>::operator=(const FunctionNode<D> &nd) {
    MWNode<D>::operator=(nd);
    return *this;
}

template<int D>
MRNode<D>* FunctionNode<D>::retrieveNode(int n, const double *r) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>* FunctionNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT;
}

/** Function evaluation.
  * Evaluate all polynomials defined on the node. */
template<int D>
double FunctionNode<D>::evalf(const double *r) {
    NOT_IMPLEMENTED_ABORT;
//    if (!this->hasCoefs()) {
//        return 0.0;
//    }
//    double arg[D];
//    int fact[D + 1];

//    MatrixXd val(this->getKp1(), D);
//    const ScalingBasis &scaling = this->tree->getScalingFunctions();

//    double scale_factor = pow(2.0, nodeIndex.scale());
//    for (int i = 0; i < D; i++) {
//        double l_factor = (double) this->nodeIndex[i];
//        arg[i] = r[i] * scale_factor - l_factor;
//    }
//    for (int i = 0; i < D + 1; i++) {
//        fact[i] = MathUtils::ipow(this->getKp1(), i);
//    }
//    scaling.evalf(arg, val);
//    double result = 0.0;
//#pragma omp parallel for shared(fact) reduction(+:result)
//    for (int i = 0; i < this->getKp1_d(); i++) {
//        double temp = (*this->coefs)(i);
//        for (int j = 0; j < D; j++) {
//            int k = (i % fact[j + 1]) / fact[j];
//            temp *= val(k, j);
//        }
//        result += temp;
//    }
//    result *= pow(2.0, 0.5 * D * nodeIndex.scale());
//    return result;
}

/*
template<int D>
void FunctionNode<D>::addCoefs(VectorXd &expCoefs, MatrixXd &mwCoefs) {
    if(not this->isAllocated()) {
    this->allocCoefs();
    }

    VectorXd &cVec = this->getCoefs();
    cVec = mwCoefs.transpose() * expCoefs;
}
*/
/** Arithmetic addition of coefficients.
  *
  * Make the function represented on THIS node the sum of the functions
  * represented on the input nodes, with constant factors. This is simply a
  * componentwise addition of coefficients.
  *
  * this->coefs = a*lhNode.coefs + b*rhNode.coefs. */
/*
template<int D>
void FunctionNode<D>::addCoefs(double a, FunctionNode<D> &lhNode,
                   double b, FunctionNode<D> &rhNode) {

    if(not this->isAllocated()) {
        this->allocCoefs();
    }

    int kp1_d = this->getKp1_d();
    int nCoefs = this->getNCoefs() - kp1_d;

    VectorXd &aVec = lhNode.getCoefs();
    VectorXd &bVec = rhNode.getCoefs();
    VectorXd &cVec = this->getCoefs();

    cVec.segment(0, kp1_d) = a*aVec.segment(0, kp1_d) + b*bVec.segment(0, kp1_d);
    VectorXd tmp = VectorXd::Zero(nCoefs);
    if (not lhNode.isGenNode()) {
    tmp += a*aVec.segment(kp1_d, nCoefs);
    }
    if (not rhNode.isGenNode()) {
    tmp += b*bVec.segment(kp1_d, nCoefs);
    }
    cVec.segment(kp1_d, nCoefs) = tmp;
}
*/
/*
template<int D>
void FunctionNode<D>::multiplyCoefs(VectorXd &expCoefs, MatrixXd &mwCoefs) {
    NOT_IMPLEMENTED_ABORT;
//    const ScalingBasis &sf = this->tree->getScalingFunctions();
//    if (sf.getType() != Interpol) {
//        NOT_IMPLEMENTED_ABORT;
//    }
//    if(not this->isAllocated()) {
//        this->allocCoefs();
//    }

//    int quadratureOrder = sf.getQuadratureOrder();
//    getQuadratureCache(qc);

//    int kp1 = this->getKp1();
//    int kp1_d = this->getKp1_d();
//    double two_scale = pow(2.0, this->getScale() + 1);

//    VectorXd sqrtWeights = qc.getWeights(quadratureOrder);
//    VectorXd sqrtWeightsInverse = sqrtWeights.array().inverse();

//    VectorXd &cVec = this->getCoefs();
//    cVec = VectorXd::Ones(this->getNCoefs());
//    double preFac = 1.0;
//    for (int n = 0; n < mwCoefs.rows(); n++) {
//        sqrtWeights = sqrtWeights.array() * sqrtWeightsInverse.array();
//        cVec = cVec.array() * mwCoefs.row(n).transpose().array();
//        preFac *= expCoefs[n];
//    }

//    cVec *= preFac;
//    sqrtWeights = two_scale * sqrtWeights.array();
//    sqrtWeights = sqrtWeights.array().sqrt();

//    int kp1_p[D];
//    for (int i = 0; i < D; i++) {
//        kp1_p[i] = MathUtils::ipow(kp1, i);
//    }

//    for (int m = 0; m < this->tDim; m++) {
//        for (int p = 0; p < D; p++) {
//            int n = 0;
//            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
//                for (int j = 0; j < kp1; j++) {
//                    for (int k = 0; k < kp1_p[p]; k++) {
//                        cVec[m * kp1_d + n] *= sqrtWeights[j];
//                        n++;
//                    }
//                }
//            }
//        }
//    }
}
*/

/** Arithmetic multiplication of coefficients.
  *
  * Make the function represented on THIS node the product of the functions
  * represented on the input nodes, with constant factors. This is NOT a
  * componentwise multiplication of coefficients, it also contains coupling of
  * quadrature weights. Only implemented for Interpolating wavelets.
  *
  * this->coefs = a*lhNode.coefs * b*rhNode.coefs. */
/*
template<int D>
void FunctionNode<D>::multiplyCoefs(double a, FunctionNode<D> &lhNode,
                    double b, FunctionNode<D> &rhNode) {
    NOT_IMPLEMENTED_ABORT;
//    const ScalingBasis &sf = this->tree->getScalingFunctions();
//    if (sf.getType() != Interpol) {
//        NOT_IMPLEMENTED_ABORT;
//    }
//    if(not this->isAllocated()) {
//        this->allocCoefs();
//    }

//    VectorXd &aVec = lhNode.getCoefs();
//    VectorXd &bVec = rhNode.getCoefs();
//    VectorXd &cVec = this->getCoefs();

//    int quadratureOrder = sf.getQuadratureOrder();
//    getQuadratureCache(qc);

//    cVec = aVec.array() * bVec.array();
//    cVec *= a*b;

//    int kp1 = this->getKp1();
//    int kp1_d = this->getKp1_d();
//    double two_scale = pow(2.0, this->getScale() + 1);

//    const VectorXd &weights = qc.getWeights(quadratureOrder);
//    VectorXd sqrtWeights = two_scale*weights.array().inverse();
//    sqrtWeights = sqrtWeights.array().sqrt();

//    int kp1_p[D];
//    for (int i = 0; i < D; i++) {
//        kp1_p[i] = MathUtils::ipow(kp1, i);
//    }

//    for (int m = 0; m < this->tDim; m++) {
//        for (int p = 0; p < D; p++) {
//            int n = 0;
//            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
//                for (int j = 0; j < kp1; j++) {
//                    for (int k = 0; k < kp1_p[p]; k++) {
//                        cVec[m * kp1_d + n] *= sqrtWeights[j];
//                        n++;
//                    }
//                }
//            }
//        }
//    }
}
*/
/** Arithmetic multiplication by a constant.
  * Componentwise multiplication of coefficients by a constant. */
template<int D>
void FunctionNode<D>::operator*=(double c) {
    if (not this->isForeign() or this->hasCoefs()) {
        this->getCoefs() *= c;
        this->calcNorms();  //maybe not the most efficient way but safest....
    }
}

/** Function integration.
  *
  * Wrapper for function integration, that requires different methods depending
  * on scaling type. Integrates the function represented on the node on the
  * full support of the node. */
template<int D>
double FunctionNode<D>::integrate() {
    if (this->isForeign()) {
        return 0.0;
    }
    switch (getFuncTree().getScalingType()) {
    case Legendre:
        return integrateLegendre();
        break;
    case Interpol:
        return integrateInterpolating();
        break;
    default:
        MSG_FATAL("Invalid scalingType");
    }
}

/** Function integration, Legendre basis.
  *
  * Integrates the function represented on the node on the full support of the
  * node. The Legendre basis is particularly easy to integrate, as the work is
  * already done when calculating its coefficients. The coefficients of the
  * node is defined as the projection integral
  *          s_i = int f(x)phi_i(x)dx
  * and since the first Legendre function is the constant 1, the first
  * coefficient is simply the integral of f(x). */
template<int D>
double FunctionNode<D>::integrateLegendre() {
    double result = this->getCoefs()[0];
    double n = (D*this->nodeIndex.getScale())/2.0;
    result *= pow(2.0, -n);
    return result;
}

/** Function integration, Interpolating basis.
  *
  * Integrates the function represented on the node on the full support of the
  * node. A bit more involved than in the Legendre basis, as is requires some
  * coupling of quadrature weights. */
template<int D>
double FunctionNode<D>::integrateInterpolating() {
    NOT_IMPLEMENTED_ABORT;
//    const ScalingBasis &sf = this->tree->getScalingFunctions();
//    VectorXd coefs = this->getCoefs();

//    int quadratureOrder = sf.getQuadratureOrder();
//    getQuadratureCache(qc);
//    const VectorXd &weights = qc.getWeights(quadratureOrder);
//    double sqWeights[this->getKp1()];
//    for (int i = 0; i < this->getKp1(); i++) {
//        sqWeights[i] = sqrt(weights[i]);
//    }

//    int kp1_p[D];
//    for (int i = 0; i < D; i++) {
//        kp1_p[i] = MathUtils::ipow(this->getKp1(), i);
//    }

//    double result = 0.0;
//    for (int p = 0; p < D; p++) {
//        int n = 0;
//        for (int i = 0; i < kp1_p[D - p - 1]; i++) {
//            for (int j = 0; j < this->getKp1(); j++) {
//                for (int k = 0; k < kp1_p[p]; k++) {
//                    coefs[n] *= sqWeights[j];
//                    n++;
//                }
//            }
//        }
//    }
//    for (int i = 0; i < this->getKp1_d(); i++) {
//        result += coefs[i];
//    }
//    double n = (D * this->nodeIndex.getScale()) / 2.0;
//    result *= pow(2.0, -n);
//    return result;
}

/** Inner product of the functions represented by the scaling basis of the nodes.
  *
  * Integrates the product of the functions represented by the scaling basis on
  * the node on the full support of the nodes. The scaling basis is fully
  * orthonormal, and the inner product is simply the dot product of the
  * coefficient vectors. Assumes the nodes have identical support. */
template<int D>
double FunctionNode<D>::scalingInnerProduct(FunctionNode<D> &inpNode) {
    int kp1_d = this->getKp1_d();
    if (this->isForeign()) {
        return 0.0;
    }
    assert(this->hasCoefs() and inpNode.hasCoefs());
    VectorXd &a = this->getCoefs();
    VectorXd &b = inpNode.getCoefs();
#ifdef HAVE_BLAS
    return cblas_ddot(kp1_d, a.data(), 1, b.data(), 1);
#else
    return a.segment(0, kp1_d).dot(b.segment(0, kp1_d));
#endif
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
  *
  * Integrates the product of the functions represented by the wavelet basis on
  * the node on the full support of the nodes. The wavelet basis is fully
  * orthonormal, and the inner product is simply the dot product of the
  * coefficient vectors. Assumes the nodes have identical support. */
template<int D>
double FunctionNode<D>::waveletInnerProduct(FunctionNode<D> &inpNode) {
    if (inpNode.isGenNode()) {
        return 0.0;
    }
    int kp1_d = this->getKp1_d();
    int nCoefs = (this->getTDim() - 1) * kp1_d;
    if (this->isForeign()) {
        return 0.0;
    }
    assert(this->hasCoefs() and inpNode.hasCoefs());
    VectorXd &a = this->getCoefs();
    VectorXd &b = inpNode.getCoefs();
#ifdef HAVE_BLAS
    return cblas_ddot(nCoefs, a.segment(kp1_d, nCoefs).data(), 1,
                      b.segment(kp1_d, nCoefs).data(), 1);
#else
    return a.segment(kp1_d, nCoefs).dot(b.segment(kp1_d, nCoefs));
#endif
}

template class FunctionNode<1>;
template class FunctionNode<2>;
template class FunctionNode<3>;

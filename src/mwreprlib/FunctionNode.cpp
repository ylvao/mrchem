#include <Eigen/Core>

#include "FunctionNode.h"
#include "FunctionTree.h"
#include "MathUtils.h"
#include "ScalingBasis.h"

#ifdef HAVE_BLAS
extern "C" {
#include BLAS_H
}
#endif

using namespace std;
using namespace Eigen;

template<int D>
FunctionNode<D>::FunctionNode(FunctionTree<D> &t, const NodeIndex<D> &nIdx)
        : MWNode<D>(t, nIdx) {
    this->squareNorm = 0.0;
}

template<int D>
FunctionNode<D>::FunctionNode(FunctionNode<D> &p, int cIdx)
        : MWNode<D>(p, cIdx) {
    this->squareNorm = 0.0;
}

template<int D>
FunctionNode<D>::FunctionNode(const MWNode<D> &n)
        : MWNode<D>(n) {
    this->squareNorm = 0.0;
}

template<int D>
FunctionNode<D>::FunctionNode(const FunctionNode<D> &n)
        : MWNode<D>(n) {
    this->squareNorm = 0.0;
}

/** Function evaluation.
  * Evaluate all polynomials defined on the node. */
template<int D>
double FunctionNode<D>::evalf(const double *r) {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");
    SET_NODE_LOCK();
    if (this->isLeafNode()) {
        this->genChildren();
        this->giveChildrenCoefs();
    }
    UNSET_NODE_LOCK();
    int cIdx = this->getChildIndex(r);
    assert(this->children[cIdx] != 0);
    return getFuncChild(cIdx).evalScaling(r);
}

template<int D>
double FunctionNode<D>::evalScaling(const double *r) const {
    if (not this->hasCoefs()) MSG_ERROR("Evaluating node without coefs");

    double arg[D];
    double n_factor = pow(2.0, this->getScale());
    const int *l_factor = this->getTranslation();
    for (int i = 0; i < D; i++) {
        arg[i] = r[i] * n_factor - (double) l_factor[i];
    }

    int fact[D + 1];
    for (int i = 0; i < D + 1; i++) {
        fact[i] = MathUtils::ipow(this->getKp1(), i);
    }

    MatrixXd val(this->getKp1(), D);
    const ScalingBasis &basis = this->getMWTree().getMRA().getScalingBasis();
    basis.evalf(arg, val);

    double result = 0.0;
#pragma omp parallel for shared(fact) reduction(+:result)
    for (int i = 0; i < this->getKp1_d(); i++) {
        double temp = (*this->coefs)(i);
        for (int j = 0; j < D; j++) {
            int k = (i % fact[j + 1]) / fact[j];
            temp *= val(k, j);
        }
        result += temp;
    }
    result *= pow(2.0, 0.5 * D * this->getScale());
    return result;
}


/** Function integration.
  *
  * Wrapper for function integration, that requires different methods depending
  * on scaling type. Integrates the function represented on the node on the
  * full support of the node. */
template<int D>
double FunctionNode<D>::integrate() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isForeign()) {
//        return 0.0;
//    }
//    if (this->isCommon() and this->tree->getRankId() != 0) {
//        return 0.0;
//    }
//    switch (getFuncTree().getScalingType()) {
//    case Legendre:
//        return integrateLegendre();
//        break;
//    case Interpol:
//        return integrateInterpolating();
//        break;
//    default:
//        MSG_FATAL("Invalid scalingType");
//    }
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
    NOT_IMPLEMENTED_ABORT;
    double result = this->getCoefs()[0];
    double n = (D*this->getScale())/2.0;
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
//    const ScalingBasis &sf = this->getMWTree().getScalingFunctions();
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
//    double n = (D * this->getScale()) / 2.0;
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
double FunctionNode<D>::dotScaling(FunctionNode<D> &inpNode) {
    NOT_IMPLEMENTED_ABORT;
//    int kp1_d = this->getKp1_d();
//    if (this->isForeign()) {
//        return 0.0;
//    }
//    assert(this->hasCoefs() and inpNode.hasCoefs());
//    VectorXd &a = this->getCoefs();
//    VectorXd &b = inpNode.getCoefs();
//#ifdef HAVE_BLAS
//    return cblas_ddot(kp1_d, a.data(), 1, b.data(), 1);
//#else
//    return a.segment(0, kp1_d).dot(b.segment(0, kp1_d));
//#endif
}

/** Inner product of the functions represented by the wavelet basis of the nodes.
  *
  * Integrates the product of the functions represented by the wavelet basis on
  * the node on the full support of the nodes. The wavelet basis is fully
  * orthonormal, and the inner product is simply the dot product of the
  * coefficient vectors. Assumes the nodes have identical support. */
template<int D>
double FunctionNode<D>::dotWavelet(FunctionNode<D> &inpNode) {
    NOT_IMPLEMENTED_ABORT;
//    if (inpNode.isGenNode()) {
//        return 0.0;
//    }
//    int kp1_d = this->getKp1_d();
//    int nCoefs = (this->getTDim() - 1) * kp1_d;
//    if (this->isForeign()) {
//        return 0.0;
//    }
//    assert(this->hasCoefs() and inpNode.hasCoefs());
//    VectorXd &a = this->getCoefs();
//    VectorXd &b = inpNode.getCoefs();
//#ifdef HAVE_BLAS
//    return cblas_ddot(nCoefs, a.segment(kp1_d, nCoefs).data(), 1,
//                      b.segment(kp1_d, nCoefs).data(), 1);
//#else
//    return a.segment(kp1_d, nCoefs).dot(b.segment(kp1_d, nCoefs));
//#endif
}

template class FunctionNode<1>;
template class FunctionNode<2>;
template class FunctionNode<3>;

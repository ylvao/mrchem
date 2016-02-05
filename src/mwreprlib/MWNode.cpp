/**
 *  Simple n-dimensional node
 *
 *  Created on: Oct 12, 2009
 *      Author: jonas
 */

#include "MWNode.h"
#include "MWTree.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

/** MWNode rootnode constructor.
  * Creates an empty rootnode given its tree and node index */
template<int D>
MWNode<D>::MWNode(MWTree<D> &t, const NodeIndex<D> &nIdx)
        : MRNode<D>(t, nIdx),
          squareNorm(-1.0),
          coefs(0) {
    clearNorms();
}

/** MWNode constructor.
  * Creates an empty node given its parent and child index */
template<int D>
MWNode<D>::MWNode(MWNode<D> &p, int cIdx)
        : MRNode<D>(p, cIdx),
          squareNorm(-1.0),
          coefs(0) {
    clearNorms();
}

template<int D>
MWNode<D>::MWNode(const MWNode<D> &n)
        : MRNode<D>(n),
          squareNorm(-1.0),
          coefs(0) {
    clearNorms();
}

/** MWNode destructor.
  * Recursive deallocation of a node and all its decendants */
template<int D>
MWNode<D>::~MWNode() {
    freeCoefs();
}

/** Allocate the coefs vector. If it is already allocated, clear the
  * HasCoefs flag and reallocate if necessary. */
template<int D>
void MWNode<D>::allocCoefs(int nCoefs) {
    if (nCoefs < 0) {
        nCoefs = this->getTDim() * this->getKp1_d();
    }
    if (this->coefs == 0) {
        this->coefs = new VectorXd(nCoefs);
    } else {
        if (this->coefs->rows() != nCoefs) { // reallocate
            delete this->coefs;
            this->coefs = new VectorXd(nCoefs);
        }
    }
    this->setIsAllocated();
    this->clearHasCoefs();
}

/** Deallocation of coefficients. */
template<int D>
void MWNode<D>::freeCoefs() {
    if (this->coefs != 0) delete this->coefs;
    this->coefs = 0;
    this->clearHasCoefs();
    this->clearIsAllocated();
}

template<int D>
void MWNode<D>::zeroCoefs() {
    if (not this->isAllocated()) {
        allocCoefs();
    }
    this->coefs->setZero();
    this->zeroNorms();
    this->setHasCoefs();
}

/** Set coefficients of node.
  *
  * Copies the argument vector to the coefficient vector of the node. Allocates
  * coefficients if needed. ASSUMES that the given vector does not exceed
  * the allocated memory of the node, and if it is smaller, trailing zeros are
  * added. Option to lock node (only used in the GenNode version of the
  * routine). */
template<int D>
void MWNode<D>::setCoefs(const Eigen::VectorXd &c) {
    if (not this->isAllocated()) {
        allocCoefs();
    }
    int nNew = c.size();
    assert (nNew <= this->getNCoefs());
    if (nNew < this->getNCoefs()) {
        this->coefs->segment(nNew, this->getNCoefs() - nNew).setZero();
    }
    this->coefs->segment(0, nNew) = c;
    this->setHasCoefs();
    this->clearNorms();
}

template<int D>
void MWNode<D>::giveChildrenCoefs(bool overwrite) {
    assert(this->isBranchNode());
    if (not this->hasCoefs()) MSG_FATAL("No coefficients!");

    ProjectedNode<D> copy(*this);
    copy.mwTransform(Reconstruction);
    const VectorXd &c = copy.getCoefs();

    int kp1_d = this->getKp1_d();
    for (int i = 0; i < this->getTDim(); i++) {
        MWNode<D> &child = this->getMWChild(i);
        if (not child.hasCoefs()) {
            child.setCoefs(c.segment(i*kp1_d, kp1_d));
        } else if (overwrite) {
            child.getCoefs().segment(0, kp1_d) = c.segment(i*kp1_d, kp1_d);
        } else {
            child.getCoefs().segment(0, kp1_d) += c.segment(i*kp1_d, kp1_d);
        }
        child.calcNorms();
    }
}

/** Takes the scaling coefficients of the children and stores them consecutively
  * in the  given vector. */
template<int D>
void MWNode<D>::copyCoefsFromChildren(VectorXd &c) {
    int kp1_d = this->getKp1_d();
    assert(this->children != 0);
    for (int i = 0; i < this->getTDim(); i++) {
        MWNode<D> &child = getMWChild(i);
        if (child.hasCoefs()) {
            VectorXd &cc = child.getCoefs();
            c.segment(i * kp1_d, kp1_d) = cc.segment(0, kp1_d);
        } else {
            c.segment(i * kp1_d, kp1_d).setZero();
        }
    }
}


/** Coefficient-Value transform
  *
  * This routine transforms the scaling coefficients of the node to the
  * function values in the corresponding quadrature roots (of its children).
  * Input parameter = forward: coef->value.
  * Input parameter = backward: value->coef.
  *
  * NOTE: this routine assumes a 0/1 (scaling on children 0 and 1)
  *       representation, in oppose to s/d (scaling and wavelet). */
template<int D>
void MWNode<D>::cvTransform(int operation) {
    NOT_IMPLEMENTED_ABORT;
//    const ScalingBasis &sf = this->getMWTree().getScalingFunctions();
//    if (sf.getType() != Interpol) {
//        NOT_IMPLEMENTED_ABORT;
//    }

//    int quadratureOrder = sf.getQuadratureOrder();
//    getQuadratureCache(qc);

//    double two_scale = pow(2.0, this->getScale() + 1);
//    VectorXd modWeights = qc.getWeights(quadratureOrder);
//    switch (operation) {
//    case Forward:
//        modWeights = modWeights.array().inverse();
//        modWeights *= two_scale;
//        modWeights = modWeights.array().sqrt();
//        break;
//    case Backward:
//        modWeights *= 1.0/two_scale;
//        modWeights = modWeights.array().sqrt();
//        break;
//    default:
//        MSG_FATAL("Invalid operation");
//    }

//    VectorXd &coefs = this->getCoefs();

//    int kp1 = this->getKp1();
//    int kp1_d = this->getKp1_d();
//    int kp1_p[D];
//    for (int d = 0; d < D; d++) {
//        kp1_p[d] = MathUtils::ipow(kp1, d);
//    }

//    for (int m = 0; m < this->getTDim(); m++) {
//        for (int p = 0; p < D; p++) {
//            int n = 0;
//            for (int i = 0; i < kp1_p[D - p - 1]; i++) {
//                for (int j = 0; j < kp1; j++) {
//                    for (int k = 0; k < kp1_p[p]; k++) {
//                        coefs[m * kp1_d + n] *= modWeights[j];
//                        n++;
//                    }
//                }
//            }
//        }
//    }
}

/** Multiwavelet transform: fast version
  *
  * Application of the filters on one node to pass from a 0/1 (scaling
  * on children 0 and 1) representation to an s/d (scaling and
  * wavelet) representation. Bit manipulation is used in order to
  * determine the correct filters and whether to apply them or just
  * pass to the next couple of indexes. The starting coefficients are
  * preserved until the application is terminated, then they are
  * overwritten. With minor modifications this code can also be used
  * for the inverse mw transform (just use the transpose filters) or
  * for the application of an operator (using A, B, C and T parts of an
  * operator instead of G1, G0, H1, H0). This is the version where the
  * three directions are operated one after the other. Although this
  * is formally faster than the other algorithm, the separation of the
  * three dimensions prevent the possibility to use the norm of the
  * operator in order to discard a priori negligible contributions.

  * Luca Frediani, August 2006
  * C++ version: Jonas Juselius, September 2009 */
template<int D>
void MWNode<D>::mwTransform(int operation) {
    int kp1 = this->getKp1();
    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
    int kp1_d = this->getKp1_d();
    const MWFilter &filter = getMWTree().getMRA().getFilter();
    VectorXd &result = getMWTree().getTmpMWCoefs();
    bool overwrite = true;

    for (int i = 0; i < D; i++) {
        int mask = 1 << i;
        for (int gt = 0; gt < this->getTDim(); gt++) {
            double *out = result.data() + gt * kp1_d;
            for (int ft = 0; ft < this->getTDim(); ft++) {
                /* Operate in direction i only if the bits along other
                 * directions are identical. The bit of the direction we
                 * operate on determines the appropriate filter/operator */
                if ((gt | mask) == (ft | mask)) {
                    double *in = this->coefs->data() + ft * kp1_d;
                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
                    MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
                    overwrite = false;
                }
            }
            overwrite = true;
        }
        this->coefs->swap(result);
    }
}

/** Set all norms to Undefined. */
template<int D>
void MWNode<D>::clearNorms() {
    this->squareNorm = -1.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = -1.0;
    }
}

/** Set all norms to zero. */
template<int D>
void MWNode<D>::zeroNorms() {
    this->squareNorm = 0.0;
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = 0.0;
    }
}

/** Calculate and store square norm and component norms, if allocated. */
template<int D>
void MWNode<D>::calcNorms() {
    this->squareNorm = calcSquareNorm();
    for (int i = 0; i < this->getTDim(); i++) {
        this->componentNorms[i] = calcComponentNorm(i);
    }
}

/** Calculate, store and return square norm. */
template<int D>
double MWNode<D>::calcSquareNorm() const {
    assert(this->isAllocated());
    assert(this->hasCoefs());
    return this->coefs->squaredNorm();
}

/** Calculate and return wavelet norm. */
template<int D>
double MWNode<D>::calcWaveletNorm() const {
    assert(this->isAllocated());
    assert(this->hasCoefs());
    int nCoefs = this->getNCoefs();
    int kp1_d = this->getKp1_d();
    return this->coefs->segment(kp1_d, nCoefs - kp1_d).norm();
}

/** Calculate the norm of one component (NOT the squared norm!). */
template<int D>
double MWNode<D>::calcComponentNorm(int i) const {
    assert(this->isAllocated());
    assert(this->hasCoefs());
    const VectorXd &c = this->getCoefs();
    int kp1_d = this->getKp1_d();
    return c.segment(i*kp1_d, kp1_d).norm();
}

template<int D>
double MWNode<D>::estimateError(bool absPrec) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isForeign()) {
//        return 0.0;
//    }
//    if (this->isCommon() and this->tree->getRankId() != 0) {
//        return 0.0;
//    }
//    double tNorm = 1.0;
//    if (not absPrec) {
//        tNorm = sqrt(getMWTree().getSquareNorm());
//    }

//    int n = this->getScale();
//    double expo = (1.0 * (n + 1));
//    double scaleFactor = max(2.0* MachinePrec, pow(2.0, -expo));
//    double wNorm = this->calcWaveletNorm();
//    double error = scaleFactor * wNorm / tNorm;
//    return error*error;
}

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. */
template<int D>
void MWNode<D>::reCompress(bool overwrite) {
    if ((not this->isGenNode()) and this->isBranchNode()) {
        if (not this->isAllocated()) {
            // This happens for seeded nodes and on distributed trees
            allocCoefs();
        }
        if (overwrite) {
            copyCoefsFromChildren(*this->coefs);
            mwTransform(Compression);
        } else {
            MatrixXd tmp = getCoefs();
            copyCoefsFromChildren(*this->coefs);
            mwTransform(Compression);
            getCoefs() += tmp;
        }
        this->setHasCoefs();
        clearNorms();
    }
}

/** Recurse down until an EndNode is found, and then crop children with
  * too high precision. */
template<int D>
bool MWNode<D>::crop(double prec, NodeIndexSet *cropIdx) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isEndNode()) {
//        return true;
//    } else {
//        assert(children != 0);
//        for (int i = 0; i < this->tDim; i++) {
//            MWNode<D> &child = *this->children[i];
//            if (child.cropChildren(prec, cropIdx)) {
//                if (not this->isForeign()) {
//                    if (this->splitCheck(prec) == false) {
//                        if (cropIdx != 0) {
//                            cropIdx->insert(&this->getNodeIndex());
//                        } else {
//                            this->deleteChildren();
//                        }
//                        return true;
//                    }
//                }
//            }
//        }
//    }
//    return false;
}

template<int D>
mpi::request MWNode<D>::isendCoefs(int who, int tag, int comp) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->hasCoefs());
//#ifdef HAVE_MPI
//    int nSend = this->getNCoefs();
//    const double *data = this->coefs->data();
//    if (comp > 0) {
//        assert(comp >= 0 and comp < this->getTDim());
//        nSend = this->getKp1_d();
//        data = data + comp * this->getKp1_d();
//    }
//    return node_group.isend(who, tag, data, nSend);
//#else
//    mpi::request dummy = 0;
//    return dummy;
//#endif
}

template<int D>
mpi::request MWNode<D>::ireceiveCoefs(int who, int tag, int comp) {
    NOT_IMPLEMENTED_ABORT;
//#ifdef HAVE_MPI
//    if (not this->isAllocated()) {
//        allocCoefs();
//    }
//    int nRecv = this->getNCoefs();
//    double *data = this->coefs->data();
//    if (comp > 0) {
//        assert(comp >= 0 and comp < this->getTDim());
//        nRecv = this->getKp1_d();
//        data = data + comp * this->getKp1_d();
//    }
//    this->setHasCoefs();
//    return node_group.irecv(who, tag, data, nRecv);
//#else
//    mpi::request dummy = 0;
//    return dummy;
//#endif
}

template class MWNode<1>;
template class MWNode<2>;
template class MWNode<3>;

/**
 *  Simple n-dimensional node
 *
 *  Created on: Oct 12, 2009
 *      Author: jonas
 */

#include "MWNode.h"
#include "MWTree.h"
#include "QuadratureCache.h"
#include "NodeIndex.h"
#include "MathUtils.h"

using namespace std;
using namespace Eigen;

/** MWNode default constructor.
  * Creates an empty node  */
template<int D>
MWNode<D>::MWNode() {
    NOT_IMPLEMENTED_ABORT;
}

/** MWNode rootnode constructor.
  * Creates an empty rootnode given its tree and node index */
template<int D>
MWNode<D>::MWNode(MRTree<D> &t, const NodeIndex<D> &nIdx) : MRNode<D>(t, nIdx) {
    this->weight[0] = 0.0;
    this->weight[1] = 0.0;
    this->componentNorms = 0;
    this->coefs = 0;
    clearNorms();
}

/** MWNode constructor.
  * Creates an empty node given its parent and child index */
template<int D>
MWNode<D>::MWNode(MWNode<D> *p, int cIdx) : MRNode<D>(p, cIdx) {
    NOT_IMPLEMENTED_ABORT;
    this->weight[0] = 0.0;
    this->weight[1] = 0.0;
    this->componentNorms = 0;
    this->coefs = 0;
    clearNorms();
}

/** MWNode copy constructor.
  *
  * Make a detached copy of a node that is not accessible through any tree. The
  * tree is still accessible from the node, as much of the node parameters are
  * in fact stored in the tree. Copying coefficients is optional. Children
  * nodes are NOT copied recursively. */
template<int D>
MWNode<D>::MWNode(const MWNode<D> &nd, MWTree<D> *t) {
    NOT_IMPLEMENTED_ABORT;
//    if (_tree == 0) {
//        MSG_FATAL("Tree not defined!");
//    }
//    parent = 0;
//    tree = _tree;
//    nodeIndex = nd.nodeIndex;
//    status = nd.status;
//    hilbertPath = nd.hilbertPath;
//    weight[0] = nd.weight[0];
//    weight[1] = nd.weight[1];
//    children = 0;
//    coefs = 0;
//    componentNorms = 0;

//    if (copyCoefs) {
//        if (nd.hasCoefs()) {
//            coefs = new VectorXd(*nd.coefs);
//            this->setHasCoefs(true);
//        }
//    } else if (nd.isAllocated()) {
//        allocCoefs(nd.getNCoefs());
//    }
//    squareNorm = nd.squareNorm;
//    if (nd.componentNorms != 0) {
//        allocComponentNorms();
//        for (int i = 0; i < tDim; i++) {
//           assert(nd.componentNorms[i] >= 0.0);
//           componentNorms[i] = nd.componentNorms[i];
//        }
//    }
//    setIsLeafNode();

//#ifdef OPENMP
//    omp_init_lock(&node_lock);
//#endif
}

/** MWNode copy constructor.
  * Make a copy of a node and assign it to another parent.
  * Copying coefficients is optional. */
template<int D>
MWNode<D>::MWNode(const MWNode<D> &nd, MWNode<D> *p) {
    NOT_IMPLEMENTED_ABORT;
//    if (_parent == 0) {
//        MSG_FATAL("Parent node not defined!");
//    }
//    parent = _parent;
//    tree = parent->tree;
//    nodeIndex = nd.nodeIndex;
//    status = nd.status;
//    hilbertPath = nd.hilbertPath;
//    weight[0] = nd.weight[0];
//    weight[1] = nd.weight[1];
//    children = 0;
//    coefs = 0;
//    componentNorms = 0;

//    if (copyCoefs) {
//        if (nd.hasCoefs()) {
//            coefs = new VectorXd(*nd.coefs);
//            this->setHasCoefs(true);
//        }
//    } else if (nd.isAllocated()) {
//        allocCoefs(nd.getNCoefs());
//    }
//    squareNorm = nd.squareNorm;
//    if (nd.componentNorms != 0) {
//        allocComponentNorms();
//        for (int i = 0; i < tDim; i++) {
//            assert(nd.componentNorms[i] >= 0.0);
//            componentNorms[i] = nd.componentNorms[i];
//        }
//    }

//    // MWNodes are abstract, so we defer the copying of children
//    if (nd.isBranchNode()) {
//        allocKindergarten();
//    }
//#ifdef OPENMP
//    omp_init_lock(&node_lock);
//#endif
}

/** MWNode equals operator.
  * Copying the content of a node, not its location.
  * Recursive copying of children is done in derived classes. */
template<int D>
MWNode<D> &MWNode<D>::operator=(const MWNode<D> &nd) {
    NOT_IMPLEMENTED_ABORT;
//    if (&nd == this) {
//        return *this;
//    }
//    if (this->tree == 0) {
//        MSG_FATAL("Cannot assign node without tree!");
//    }
//    if (not this->checkCompatible(nd)) {
//        MSG_FATAL("Nodes not compatible!");
//    }

//    SET_NODE_LOCK();

//    // MWNodes are abstract, so we defer the copying of children
//    // NB! This must be done _before_ copying, not to mess up the status flags
//    if (this->isBranchNode()) {
//        deleteChildren();
//    }
//    status = nd.status;
//    if (nd.isBranchNode()) {
//        allocKindergarten();
//    }
//    nodeIndex = nd.nodeIndex;
//    hilbertPath = nd.hilbertPath;
//    squareNorm = nd.squareNorm;
//    weight[0] = nd.weight[0];
//    weight[1] = nd.weight[1];
//    if (nd.componentNorms != 0) {
//        allocComponentNorms();
//        for (int i = 0; i < tDim; i++) {
//            assert(nd.componentNorms[i] != Undefined);
//            componentNorms[i] = nd.componentNorms[i];
//        }
//    }

//    MWNode<D>::clearCoefs(); // clear any existing coefs
//    if (nd.isAllocated()) {
//        if (nd.hasCoefs()) {
//            this->coefs = new VectorXd(*nd.coefs);
//            this->setHasCoefs(true);
//        } else {
//            this->coefs = new VectorXd(nd.getNCoefs());
//            this->setIsAllocated();
//        }
//    }
//    UNSET_NODE_LOCK();

//    return *this;
}

/** MWNode destructor.
  * Recursive deallocation of a node and all its decendants */
template<int D>
MWNode<D>::~MWNode() {
    freeCoefs();
    freeComponentNorms();
}

/** Recurse down until an EndNode is found, and then crop children with
  * too high precision. */
template<int D>
bool MWNode<D>::cropChildren(double prec, set<const NodeIndex<D> *,
                             NodeIndexComp<D> > *cropIdx) {
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
    this->setHasCoefs(false);
}

/** Allocate the coefs vector with nCoefs coefficients. If it is already
  * allocated, keep the overlapping coefficients, set any additional new coefs
  * to zero or delete any extra old coefs. */
template<int D>
void MWNode<D>::reallocCoefs(int nCoefs) {
    NOT_IMPLEMENTED_ABORT;
    if (nCoefs != this->getNCoefs()) {
        VectorXd *oldCoefs = this->coefs;
        this->coefs = 0;
        allocCoefs(nCoefs);
        if (this->getNCoefs() >= oldCoefs->size()) {
            setCoefs(*oldCoefs);
        } else {
            setCoefs(oldCoefs->segment(0, this->getNCoefs()));
        }
        delete oldCoefs;
    }
}

/** Deallocation of coefficients. */
template<int D>
void MWNode<D>::freeCoefs() {
    if (this->coefs != 0) {
        delete this->coefs;
    }
    this->setHasCoefs(false);
    this->clearIsAllocated();
    this->coefs = 0;
}

template<int D>
void MWNode<D>::zeroCoefs() {
    assert(this->coefs != 0);
    this->coefs->setZero();
    this->setHasCoefs(true);
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

//    double two_scale = pow(2.0, getScale() + 1);
//    VectorXd modWeights = qc.getWeights(quadratureOrder);
//    switch (operation) {
//        case Forward:
//            modWeights = modWeights.array().inverse();
//            modWeights *= two_scale;
//            modWeights = modWeights.array().sqrt();
//            break;
//        case Backward:
//            modWeights *= 1.0/two_scale;
//            modWeights = modWeights.array().sqrt();
//            break;
//        default:
//            MSG_FATAL("Invalid operation");
//    }

//    VectorXd &coefs = getCoefs();

//    int kp1 = getKp1();
//    int kp1_d = getKp1_d();
//    int kp1_p[D];
//    for (int i = 0; i < D; i++) {
//        kp1_p[i] = MathUtils::ipow(getKp1(), i);
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
    NOT_IMPLEMENTED_ABORT;
//    int kp1 = this->getKp1();
//    int kp1_dm1 = MathUtils::ipow(kp1, D - 1);
//    int kp1_d = this->getKp1_d();
//    const Filter &filter = getMWTree().getFilter();
//    VectorXd &result = getMWTree().getTmpMWCoefs();
//    bool overwrite = true;

//    for (int i = 0; i < D; i++) {
//        int mask = 1 << i;
//        for (int gt = 0; gt < this->getTDim(); gt++) {
//            double *out = result.data() + gt * kp1_d;
//            for (int ft = 0; ft < this->getTDim(); ft++) {
//                /* Operate in direction i only if the bits along other
//                 * directions are identical. The bit of the direction we
//                 * operate on determines the appropriate filter/operator */
//                if ((gt | mask) == (ft | mask)) {
//                    double *in = this->coefs->data() + ft * kp1_d;
//                    int fIdx = 2 * ((gt >> i) & 1) + ((ft >> i) & 1);
//                    const MatrixXd &oper = filter.getSubFilter(fIdx, operation);
//                    MathUtils::applyFilter(out, in, oper, kp1, kp1_dm1, overwrite);
//                    overwrite = false;
//                }
//            }
//            overwrite = true;
//        }
//        this->coefs->swap(result);
//    }
}

/** Set all norms to Undefined. */
template<int D>
void MWNode<D>::clearNorms() {
    this->squareNorm = -1.0;
    if (this->componentNorms != 0) {
        for (int i = 0; i < this->getTDim(); i++) {
            this->componentNorms[i] = -1.0;
        }
    }
}

/** Deallocate component norms. */
template<int D>
void MWNode<D>::freeComponentNorms() {
    if (this->componentNorms != 0) {
        delete [] this->componentNorms;
        this->componentNorms = 0;
    }
}

/** Set all norms to zero. */
template<int D>
void MWNode<D>::zeroNorms() {
    this->squareNorm = 0.0;
    if (this->componentNorms != 0) {
        for (int i = 0; i < this->getTDim(); i++) {
            this->componentNorms[i] = 0.0;
        }
    }
}

/** Calculate and store square norm and component norms, if allocated. */
template<int D>
void MWNode<D>::calcNorms() {
    NOT_IMPLEMENTED_ABORT;
    calcSquareNorm();
    if (this->componentNorms != 0) {
        calcComponentNorms();
    }
}

/** Calculate, store and return square norm. */
template<int D>
double MWNode<D>::calcSquareNorm() {
    NOT_IMPLEMENTED_ABORT;
    assert(this->coefs != 0);
    assert(this->hasCoefs());
    this->squareNorm = this->coefs->squaredNorm();
    return this->squareNorm;
}

/** Calculate and return scaling norm. */
template<int D>
double MWNode<D>::calcScalingNorm() {
    NOT_IMPLEMENTED_ABORT;
    assert(this->hasCoefs());
    return this->coefs->segment(0, this->getKp1_d()).squaredNorm();
}

/** Calculate and return wavelet norm. */
template<int D>
double MWNode<D>::calcWaveletNorm() {
    NOT_IMPLEMENTED_ABORT;
    assert(this->hasCoefs());
    int nCoefs = this->getNCoefs();
    int kp1_d = this->getKp1_d();
    return this->coefs->segment(kp1_d, nCoefs - kp1_d).norm();
}

// Don't bother trying to inline the following functions; they are virtual

/** Return coefficients, allocate if necessary. */
template<int D>
VectorXd& MWNode<D>::getCoefs() {
    NOT_IMPLEMENTED_ABORT;
    if (not this->isAllocated()) { // Lazy allocation of the fly
        allocCoefs();
    }
    return *this->coefs;
}

/** Return coefficients, allocate if necessary. Const version.*/
template<int D>
const VectorXd& MWNode<D>::getCoefs() const {
    NOT_IMPLEMENTED_ABORT;
    assert(this->coefs != 0);
    return *this->coefs;
}

/** Takes the scaling coefficients of the children and stores them consecutively
  * in the  given vector. */
template<int D>
void MWNode<D>::copyScalingCoefsFromChildren(VectorXd &scaling) {
    NOT_IMPLEMENTED_ABORT;
//    int kp1_d = this->getKp1_d();
//    assert(this->children != 0);
//    for (int i = 0; i < this->getTDim(); i++) {
//        MWNode<D> &child = getMWChild(i);
//        if (child.hasCoefs()) {
//            VectorXd &cc = child.getCoefs();
//            scaling.segment(i * kp1_d, kp1_d) = cc.segment(0, kp1_d);
//        } else {
//            scaling.segment(i * kp1_d, kp1_d).setZero();
//        }
//    }
}

/** Update the coefficients of the node by a mw transform of the scaling
  * coefficients of the children. Option to overwrite or add up existing
  * coefficients. */
template<int D>
void MWNode<D>::reCompress(bool overwrite) {
    NOT_IMPLEMENTED_ABORT;
//    if ((not this->isGenNode()) and this->isBranchNode()) {
//        if (not this->isAllocated()) {
//            // This happens for seeded nodes and on distributed trees
//            allocCoefs();
//        }
//        if (overwrite) {
//            copyScalingCoefsFromChildren(*this->coefs);
//            mwTransform(Compression);
//        } else {
//            MatrixXd tmp = getCoefs();
//            copyScalingCoefsFromChildren(*this->coefs);
//            mwTransform(Compression);
//            getCoefs() += tmp;
//        }
//        this->setHasCoefs();
//        clearNorms();
//    }
}

/** Testing if the branch decending from this node differs from the branch
  * decending from another node. Tests coefficients, norms and status bits.
  * Returns true if branches differ. */
template<int D>
bool MWNode<D>::diffBranch(const MWNode<D> &rhs) const {
    NOT_IMPLEMENTED_ABORT;
//    bool differ = false;
//    if (this->getRankId() != rhs.getRankId()) {
//        differ = true;
//        println(1, "Nodes have different rankId!");
//        goto check_branch;
//    }
//    if (this->isForeign()) {
//        goto check_branch;
//    }
//    if (this->diffCoefs(rhs)) {
//        println(1, "Coefs differ in node: " << this->getNodeIndex());
//        differ = true;
//    }
//    if (this->diffNorms(rhs)) {
//        println(1, "Norms differ in node: " << this->getNodeIndex());
//        differ = true;
//    }
//check_branch:
//    if (this->status != rhs.status) {
//        println(1, "Status differ in node: " << this->getNodeIndex());
//        for (int i=7; i>=0; i--) {
//            int bit = ((this->status >> i) & 1);
//            printout(2, bit);
//        }
//        cout << endl;
//        for (int i=7; i>=0; i--) {
//            int bit = ((rhs.status >> i) & 1);
//            printout(1, bit);
//        }
//        println(2, "");
//        println(2, "");
//        differ = true;
//    }
//    if (differ) {
//        println(1, *this);
//        println(1, rhs);
//        println(1, "");
//    }
//    if (this->isBranchNode()) {
//        if (rhs.isBranchNode()) {
//            for (int i = 0; i < this->getTDim(); i++) {
//                if (this->children[i]->diffBranch(*rhs.children[i])) {
//                    return true;
//                }
//           }
//        } else{
//            println(2, "Children differ in node: " << this->getNodeIndex());
//            return true;
//        }
//    } else {
//        if (rhs.isBranchNode()) {
//            println(2, "Children differ in node: " << this->getNodeIndex());
//            return true;
//        }
//    }
//    return differ;
}

/** Testing if the coefficients of this node differs from the coefficients
  * of another node. Returns true if coefs differ. */
template<int D>
bool MWNode<D>::diffCoefs(const MWNode<D> &rhs) const {
    NOT_IMPLEMENTED_ABORT;
    bool differ = false;
    if (this->hasCoefs() != rhs.hasCoefs()) {
        println(2, "hasCoefs mismatch: ");
        return true;
    }
    if (this->hasCoefs() and rhs.hasCoefs()) {
        const VectorXd &lhsCoefs = this->getCoefs();
        const VectorXd &rhsCoefs = rhs.getCoefs();
        if (lhsCoefs.size() != rhsCoefs.size()) {
            println(2, "nCoefs mismatch: " << lhsCoefs.size() <<
                    "!=" << rhsCoefs.size());
            return true;
        }
        differ = false;
        for (int i=0; i < this->getNCoefs(); i++) {
            if (fabs(lhsCoefs[i] - rhsCoefs[i]) > MachinePrec) {
                differ = true;
                println(3, "Coefs differ at " << this->nodeIndex <<
                ": "<< lhsCoefs[i] << " != " << rhsCoefs[i]);
            }
        }
    }
    return differ;
}

/** Testing if the norms of this node differs from the coefficients of
  * another node. Returns true if norms differ. */
template<int D>
bool MWNode<D>::diffNorms(const MWNode<D> &rhs) const {
    NOT_IMPLEMENTED_ABORT;
    if (this->squareNorm != rhs.squareNorm) {
        println(2, "SquareNorm mismatch: ");
        println(2, this->squareNorm);
        println(2, rhs.squareNorm << endl);
        return true;
    }
    return false;
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
    NOT_IMPLEMENTED_ABORT;
    if (not this->isAllocated()) {
        allocCoefs();
    }
    int nNew = c.size();
    assert (nNew <= this->getNCoefs());
    if (nNew < this->getNCoefs()) {
        this->coefs->segment(nNew, this->getNCoefs() - nNew).setZero();
    }
    this->coefs->segment(0, nNew) = c;
    this->setHasCoefs(true);
}

/** Check if node should be seeded.
  *
  * Performing the initial seedCheck regarding max depth and uniform depth
  * of the tree, before the proper "func" dependent seedCheck is performed. */
//template<int D>
//bool MWNode<D>::seedCheck(RepresentableObject<D> &func) {
//    int scale = this->getScale();
//    int depth = this->getDepth();
//    if (scale >= this->getMWTree().getMaxScale()) {
//        println(1, "+++ Maximum depth reached: " << this->getDepth());
//        return false;
//    } else if (depth < this->getMWTree().getUniform()) {
//        return true;
//    } else {
//        return func.checkSeedNode(*this);
//    }
//}

/** Get the quadrature root matrix of all children, stored consecutively in
  * a vector. */
template<int D>
void MWNode<D>::getChildrenQuadRoots(vector<MatrixXd *> &quadPts) {
    NOT_IMPLEMENTED_ABORT;
    int kp1 = this->getKp1();
    double sfac = pow(2.0, this->getScale() + 1);

    getQuadratureCache(qc);
    const VectorXd &pts = qc.getRoots(kp1);

    for (int child = 0; child < this->getTDim(); child++) {
        MatrixXd *tmpMat = new MatrixXd(kp1, D);
        int l[D];
        this->calcChildTranslation(child, l);
        for (int d = 0; d < D; d++) {
            for (int i = 0; i < kp1; i++) {
                (*tmpMat)(i, d) = (pts(i) + l[d]) / sfac;
            }
        }
        quadPts.push_back(tmpMat);
    }
}


template class MWNode<1>;
template class MWNode<2>;
template class MWNode<3>;

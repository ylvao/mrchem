/*
 *
 *
 *  \date Oct 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif
 */

#include "GenNode.h"

using namespace std;
using namespace Eigen;

/** GenNode constructor.
  * Creates an empty node given its parent and translation vector */
template<int D>
GenNode<D>::GenNode(FunctionNode<D> &p, int cIdx) : FunctionNode<D> (p, cIdx) {
    NOT_IMPLEMENTED_ABORT;
//    this->setIsGenNode();
//    if (this->parent->isGenNode()) {
//        this->genRootNode = static_cast<GenNode<D> *>(parent)->getGenRootNode();
//    } else {
//        this->genRootNode = static_cast<ProjectedNode<D> *>(parent);
//    }
//    this->tree->incrementGenNodeCount();
}

/** GenNode equals operator.
  * Copying the content of a node, not its location. Includes recursive copying
  * of children nodes. */
template<int D>
GenNode<D> &GenNode<D>::operator=(const GenNode<D> &nd) {
    NOT_IMPLEMENTED_ABORT;
//    if (this == &node) {
//        return *this;
//    }
//    FunctionNode<D>::operator=(node);
//    if (this->tree == 0) {
//        MSG_FATAL("Cannot assign node without tree!");
//    }
//    this->tree->incrementGenNodeCount();
//    this->genRootNode = node.genRootNode;
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        for (int i = 0; i < this->tDim; i++) {
//            if (this->hasChild(i)) {
//                const GenNode<D> &child =
//                        static_cast<const GenNode<D> &> (node.getMWChild(i));
//                this->children[i] = new GenNode<D> (child, this, true);
//            }
//        }
//    }
//    return *this;
}

/** GenNode destructor.
  * Base class destructor deletes children nodes recursively. */
template<int D>
GenNode<D>::~GenNode() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->tree != 0) {
//        this->tree->decrementGenNodeCount();
//        if (this->isAllocated()) {
//            this->tree->decrementAllocGenNodeCount();
//        }
//    }
}

/** Deallocates all generated coefs recursively, but retains the norms */
template<int D>
void GenNode<D>::clearCoefs() {
    NOT_IMPLEMENTED_ABORT;
//    releaseCoefs();
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        for (int i=0; i < this->tDim; i++) {
//            if (this->children[i] != 0) {
//                this->children[i]->clearCoefs();
//            }
//        }
//    }
}

/** Allocate the coefs vector. GenNodes will by default only allocate memory
  * for the (k+1)^d scaling coefficients, and are thus 2^(-D) smaller than real
  * nodes. If it is already allocated, clear the HasCoefs flag and reallocate
  * if necessary. */
template<int D>
void GenNode<D>::allocCoefs(int nCoefs) {
    NOT_IMPLEMENTED_ABORT;
//    if (nCoefs < 0) {
//        nCoefs = this->getKp1_d();
//    }
//    if (this->coefs != 0) {
//        if (this->coefs->rows() != nCoefs) { // reallocate
//            delete this->coefs;
//            this->coefs = new VectorXd(nCoefs);
//        }
//    } else {
//        this->coefs = new VectorXd(nCoefs);
//    }
//    this->coefs->setZero();
//    this->setIsAllocated();
//    this->clearHasCoefs();
}


/** Allocating children nodes.
  *
  * This routine creates 2^D empty GenNode children nodes with the
  * appropriate translation and Hilbert path parameters. */
template<int D>
void GenNode<D>::createChildren() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->children == 0) {
//        this->allocKindergarten();
//    }
//    for (int i = 0; i < this->tDim; i++) {
//        createChild(i);
//    }
//    this->setIsBranchNode();
}

/** Allocating child node.
  *
  * Given a child index, this routine creates an empty GenNodechild node with
  * the appropriate translation and Hilbert path parameters. */
template<int D>
void GenNode<D>::createChild(int i) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->children != 0);
//    assert(this->children[i] == 0);
//    int l[D];
//    this->calcChildTranslation(i, l);
//    GenNode<D> *child = new GenNode<D> (this, l);
//    this->children[i] = child;
//    return *child;
}

/** Generating children nodes.
  *
  * This routine creates 2^D children GenNodes with the appropriate translation
  * and Hilbert path parameters. Option to calculate the children scaling coefs
  * by MW transform. */
template<int D>
void GenNode<D>::genChildren(bool genEmpty) {
    NOT_IMPLEMENTED_ABORT;
//    createChildren();
//    if (not genEmpty) {
//        if (this->hasCoefs()) {
//            this->giveChildrenScaling();
//        }
//        if (this->tree->getAutoClean()) {
//            this->releaseCoefs();
//        }
//    }
}

/** Given the scaling and wavelet coefficients of the node, do a wavelet
  * transform and distribute scaling coefficient to all children. Option to
  * overwrite or add up existing coefficients. */
template<int D>
void GenNode<D>::giveChildrenScaling(bool overwrite) {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->isBranchNode());
//    if (not this->hasCoefs()) {
//        MSG_FATAL("No coefficients! You need to project first!")
//    }
//    int kp1_d = this->getKp1_d();
//    int nCoefs = this->tDim * kp1_d;
//    GenNode<D> tmpNode(*this, this->tree, false);
//    tmpNode.allocCoefs(nCoefs);
//    tmpNode.setCoefsNoLock(*this->coefs);
//    tmpNode.mwTransform(Reconstruction);
//    for (int i = 0; i < this->tDim; i++) {
//        MWNode<D> &child = this->getMWChild(i);
//        if (not child.hasCoefs()) {
//            child.setCoefs(tmpNode.getCoefs().segment(i * kp1_d, kp1_d));
//        } else if (overwrite) {
//            child.getCoefs().segment(0, kp1_d) =
//                tmpNode.getCoefs().segment(i * kp1_d, kp1_d);
//        } else {
//            child.getCoefs().segment(0, kp1_d) +=
//                tmpNode.getCoefs().segment(i * kp1_d, kp1_d);
//        }
//        child.calcNorms();
//    }
//    tmpNode.clearTreePointer();
}

/** Distributing scaling coefs to the nodes that share the same parent.
  *
  * Assuming this node contains the scaling coefs of all 2^D siblings,
  * this routine distributes them to their rightful owners, leaving wavelet
  * coefs zero. Siblings that already have coefs are overwritten. *this GenNode
  * is already locked at this point, so one has to use the setCoefs() instead of
  * setLockedCoefs(). */
template<int D>
void GenNode<D>::giveSiblingsScaling(VectorXd &coefs) {
    NOT_IMPLEMENTED_ABORT;
//    int kp1_d = this->getKp1_d();
//    int tDim = this->tDim;
//    int myIdx = this->parent->getChildIndex(this->getNodeIndex());
//    for (int i = 0; i < tDim; i++) {
//        // If this is a GenNode, so are the siblings!
//        GenNode<D> &sibling =
//                static_cast<GenNode<D> &>(this->parent->getMWChild(i));
//        if (i != myIdx) {
//            sibling.setCoefsNoLock(coefs.segment(i*kp1_d, kp1_d));
//        } else {
//            sibling.setCoefsNoLock(coefs.segment(i*kp1_d, kp1_d));
//        }
//    }
}

/** If an existing, deallocated/purged GenNode is asked for it's coefs,
  * this routine regenerates the coefs on the fly. */
template<int D>
void GenNode<D>::regenerateCoefs() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->hasCoefs()) {
//        return;
//    }
//    int idx[this->tree->getMaxDepth()];
//    int kp1_d = this->getKp1_d();

//    ProjectedNode<D> tmpNode(*this->getGenRootNode(), &this->getMWTree());
//    int relDepth = this->getDepth() - tmpNode.getDepth();
//    makeIndexList(idx, relDepth, this->nodeIndex);
//    for (int i = 0; i < relDepth; i++) {
//        tmpNode.mwTransform(Reconstruction);
//        if (i == (relDepth - 1)) {
//            giveSiblingsScaling(tmpNode.getCoefs());
//        } else {
//            tmpNode.setCoefs(tmpNode.getCoefs().segment(idx[i] * kp1_d, kp1_d));
//        }
//    }
//    tmpNode.clearTreePointer();
}

/** Recursively make a list of the indexes needed to traverse
  * back to the current node from a node relDepth above. The routine
  * stops when it hits a proper node, or depth is zero. */
template<int D>
void GenNode<D>::makeIndexList(int *idxList, int depth, const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT;
//    depth--;
//    idxList[depth] = this->parent->getChildIndex(idx);
//    if (depth == 0) {
//        return;
//    }
//    if (this->parent->isGenNode()) {
//        static_cast<GenNode<D> *>
//            (this->parent)->makeIndexList(idxList, depth, idx);
//    }
//    return;
}

/** Routine to find a particular node. Keep recursing until the last
  * node is found. Assumes that the requested coordinate is in fact defined on
  * the first node.*/
template<int D>
MWNode<D> *GenNode<D>::retrieveNode(int n, const double *r) {
    NOT_IMPLEMENTED_ABORT;
//    lockSiblings();
//    if (not this->hasCoefs()) {
//        regenerateCoefs();
//    }
//    unlockSiblings();
//    return this;
}

/** Routine to find a particular node. Keep recursing until the node
  * is found or scale n is reached. Assumes that the requested node is in fact a
  * decendant of the first node.*/
template<int D>
MWNode<D> *GenNode<D>::retrieveNode(const NodeIndex<D> &idx, bool genEmpty) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nodeIndex.scale() == idx.getScale()) { // we're done
//        if (genEmpty) {
//            return this;
//        }
//        lockSiblings();
//        if (not this->hasCoefs()) {
//            regenerateCoefs();
//        }
//        unlockSiblings();
//        return this;
//    }
//    SET_NODE_LOCK();
//    if (this->isLeafNode()) {
//        genChildren(genEmpty);
//    }
//    UNSET_NODE_LOCK();
//    int childIdx = MWNode<D>::getChildIndex(idx);
//    assert(this->children != 0);
//    GenNode<D> *child = static_cast<GenNode<D> *>(this->children[childIdx]);
//    return child->retrieveNode(idx, genEmpty);
}

/** Deallocate coefficient vector of GenNode, after calculating the square norm
  * of the node. */
template<int D>
void GenNode<D>::releaseCoefs() {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->isAllocated()) {
//        return;
//    }
//    if (this->hasCoefs()) {
//        // Calculate norms if they have not been calculated yet
//        if (this->squareNorm < 0.0) {
//            this->calcSquareNorm();
//        }
//    }
//    MWNode<D>::clearCoefs();
//    this->tree->decrementAllocGenNodeCount();
}

/** Set coefficient vector of a GenNode.
  *
  * Copies the argument vector to the coefficient vector of the node. Allocates
  * coefficients if needed. ASSUMES that the given vector does not exceed
  * the allocated memory of the node, and if it is smaller, trailing zeros are
  * added.  */
template<int D>
void GenNode<D>::setCoefs(const Eigen::VectorXd &c) {
    NOT_IMPLEMENTED_ABORT;
//    SET_NODE_LOCK();
//    if (not this->isAllocated()) {
//        this->getMWTree().incrementAllocGenNodeCount();
//    }
//    MWNode<D>::setCoefs(c);
//    UNSET_NODE_LOCK();
}
/** Set coefficient vector of a GenNode. */
template<int D>
void GenNode<D>::setCoefsNoLock(const Eigen::VectorXd &c) {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->isAllocated()) {
//        this->getMWTree().incrementAllocGenNodeCount();
//    }
//    MWNode<D>::setCoefs(c);
}

// Don't bother trying to inline the following functions; they are virtual*/

/** Get coefficients of GenNode, regenerate if needed. */
template<int D>
VectorXd& GenNode<D>::getCoefs() {
    NOT_IMPLEMENTED_ABORT;
//    lockSiblings();
//    if (not this->hasCoefs()) {
//        regenerateCoefs();
//    }
//    unlockSiblings();
//    return MWNode<D>::getCoefs();
}

/** Get coefficients of GenNode, regenerate if needed, without locking. */
template<int D>
VectorXd& GenNode<D>::getCoefsNoLock() {
    NOT_IMPLEMENTED_ABORT;
//    if (not this->hasCoefs()) {
//        regenerateCoefs();
//    }
//    return MWNode<D>::getCoefs();
}

/** Get coefficients of GenNode, const version that does not regenerate. */
template<int D>
const VectorXd& GenNode<D>::getCoefs() const {
    NOT_IMPLEMENTED_ABORT;
//    assert(this->hasCoefs());
//    return MWNode<D>::getCoefs();
}

/** Copies the coefs from the first ancestor that HAS coefs.
  *
  * Virtually resolved to stop recursion when a ProjectedNode
  * is encountered. Returns the depth of the node it copied from. */
template<int D>
int GenNode<D>::getGenRootCoefs(VectorXd &coefs) {
    NOT_IMPLEMENTED_ABORT;
//    bool foundCoefs = false;
//    SET_NODE_LOCK();
//    if (this->hasCoefs()) {
//        coefs = VectorXd::Zero(this->tDim * this->getKp1_d());
//        coefs.segment(0, this->getNCoefs()) = *this->coefs;
//        foundCoefs = true;
//    }
//    UNSET_NODE_LOCK();

//    int depth;
//    if (foundCoefs) {
//        depth = this->getDepth();
//    } else {
//        depth = this->getFuncParent().getGenRootCoefs(coefs);
//    }
//    return depth;
}

/** Clear coefficients of generated nodes.
  *
  * The node structure is kept, only the coefficients are cleared. */
template<int D>
void GenNode<D>::clearGenerated() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        for (int i = 0; i < this->tDim; i++) {
//            if (this->children[i] != 0) {
//                this->getFuncChild(i).clearGenerated();
//            }
//        }
//    }
//    clearCoefs();
}

/** Remove all generated nodes recursively.
  *
  * All generated nodes are removed (nothing is kept). */
template<int D>
void GenNode<D>::purgeGenerated() {
    NOT_IMPLEMENTED_ABORT;
//    if (this->isBranchNode()) {
//        assert(this->children != 0);
//        this->deleteChildren();
//    }
}

/** MW transform for GenNodes.
  *
  * IMPORTANT:
  * This routine increases nCoefs of the node from kp1_d to 2^d * kp1_d.
  * This means that YOU are responsible for the deallocation. However, it should
  * not cause any problems if you forget, it will just use more memory (that is,
  * the code should support both vector sizes in the same tree).
  * If you are doing a compression, make sure your coef vector contains all
  * scaling coefs of the children, and thus is full size, otherwise the missing
  * coefs will be set to zero, and your transform will return cruft. */
template<int D>
void GenNode<D>::mwTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
//    int totCoefs = this->tDim * this->getKp1_d();
//    if (this->getNCoefs() != totCoefs) {
//        this->reallocCoefs(totCoefs);
//    }
//    MWNode<D>::mwTransform(kind);
}

/** Coefficient-value transform for GenNodes. Not implemented. */
template<int D>
void GenNode<D>::cvTransform(int kind) {
    NOT_IMPLEMENTED_ABORT;
//    int nCoefs = this->getCoefs().size();
//    int totCoefs = this->tDim * this->getKp1_d();
//    if (nCoefs != totCoefs) {
//        MSG_FATAL("Cannot do cvTransform without wavelet coefs");
//    }
//    MWNode<D>::cvTransform(kind);
}

/** Activate the OMP node lock on all 2^D sibling nodes. */
template<int D>
void GenNode<D>::lockSiblings() {
    NOT_IMPLEMENTED_ABORT;
//    MWNode<D> *parent = &this->getMWParent();
//    if (parent != 0) {
//        /* Since all threads set the locks in the same order starting from 0,
//        there is no risk of a deadlock here. */
//        for (int i = 0; i < parent->getNChildren(); i++) {
//            parent->getMWChild(i).lockNode();
//        }
//    }
}

/** Deactivate the OMP node lock on all 2^D sibling nodes. */
template<int D>
void GenNode<D>::unlockSiblings() {
    NOT_IMPLEMENTED_ABORT;
//    MWNode<D> *parent = &this->getMWParent();
//    if (parent != 0) {
//        for (int i = 0; i < parent->getNChildren(); i++) {
//            parent->getMWChild(i).unlockNode();
//        }
//    }
}

template<int D>
double GenNode<D>::evalf(const double *r) {
    NOT_IMPLEMENTED_ABORT;
//    lockSiblings();
//    if (not this->hasCoefs()) {
//        regenerateCoefs();
//    }
//    unlockSiblings();
//    return FunctionNode<D>::evalf(r);
}

template<int D>
void GenNode<D>::incrementNodeWeight(int i, double w) {
    NOT_IMPLEMENTED_ABORT;
//    this->genRootNode->incrementNodeWeight(i, w);
}

template<int D>
double GenNode<D>::getComponentNorm(int i) {
    NOT_IMPLEMENTED_ABORT;
//    assert(i >= 0 and i < this->tDim);
//    if (i != 0) {
//        return 0.0;
//    }
//    lockSiblings();
//    if (this->squareNorm < 0.0) {
//        // Since it's expensive to lock, we calculate all siblings in one go
//        MWNode<D> *parent = &this->getMWParent();
//        if (parent != 0) {
//            for (int i = 0; i < parent->getNChildren(); i++) {
//                GenNode<D> &node =
//                        static_cast<GenNode<D> &>(parent->getMWChild(i));
//                node.squareNorm = node.getCoefsNoLock().squaredNorm();
//            }
//        }
//    }
//    unlockSiblings();
//    return this->squareNorm;
}

template<int D>
double GenNode<D>::calcSquareNorm() {
    NOT_IMPLEMENTED_ABORT;
//    this->squareNorm = -1.0;
//    return this->getComponentNorm(0);
}

template<int D>
double GenNode<D>::calcScalingNorm() {
    NOT_IMPLEMENTED_ABORT;
//    this->squareNorm = -1.0;
//    return this->getComponentNorm(0);
}

template class GenNode<1> ;
template class GenNode<2> ;
template class GenNode<3> ;

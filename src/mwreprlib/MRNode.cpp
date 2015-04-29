#include "MRNode.h"
#include "MathUtils.h"

using namespace std;

template<int D>
MRNode<D>::MRNode() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>::MRNode(MRTree<D> &t, const NodeIndex<D> &nIdx) : nodeIndex(nIdx) {
    this->tree = &t;
    this->tree->incrementNodeCount(getScale());

    this->parent = 0;
    this->status = 0;
    this->children = 0;

    setIsLeafNode();
    setIsRootNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
MRNode<D>::MRNode(MRNode<D> *p, int cIdx) {
    this->parent = p;
    this->status = 0;
    this->children = 0;

    if (this->parent == 0) {
        NOT_IMPLEMENTED_ABORT;
    } else {
        this->parent->calcChildNodeIndex(this->nodeIndex, cIdx);
        this->tree = this->parent->tree;
        this->tree->incrementNodeCount(getScale());
    }
    this->setRankId(this->parent->getRankId());
    setIsLeafNode();
#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
MRNode<D>::MRNode(const MRNode<D> &nd, MRNode<D> *p) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>::MRNode(const MRNode<D> &nd, MRTree<D> *t) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>& MRNode<D>::operator=(const MRNode<D> &nd) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>::~MRNode() {
    lockNode();
    if (this->isBranchNode()) {
        assert(this->children != 0);
        deleteChildren();
    }
    this->tree->decrementNodeCount(getScale());
    unlockNode();
#ifdef OPENMP
    omp_destroy_lock(&node_lock);
#endif
}

template<int D>
void MRNode<D>::allocKindergarten() {
    assert(this->children == 0);
    int tDim = getTDim();
    this->children = new MRNode<D> *[tDim];
    for (int i = 0; i < tDim; i++) {
        this->children[i] = 0;
    }
}

/** Recurcive deallocation of children and all their decendants.
  * Leaves node as LeafNode and children[] as null pointer. */
template<int D>
void MRNode<D>::deleteChildren() {
    assert(this->children != 0);
    for (int n = 0; n < getTDim(); n++) {
        assert(this->children[n] != 0);
        delete this->children[n];
        this->children[n] = 0;
    }
    delete [] this->children;
    this->children = 0;
    this->setIsLeafNode();
}


template<int D>
void MRNode<D>::createChildren() {
    assert(this->children == 0);
    this->allocKindergarten();
    for (int n = 0; n < getTDim(); n++) {
        createChild(n);
    }
    this->setIsBranchNode();
    this->clearIsEndNode();
}

template<int D>
void MRNode<D>::genChildren() {
    assert(this->children == 0);
    this->allocKindergarten();
    int nChildren = this->getTDim();
    for (int i = 0; i < nChildren; i++) {
        genChild(i);
    }
    this->setIsBranchNode();
}

template<int D>
void MRNode<D>::purgeGenerated() {
    if (this->isBranchNode()) {
        assert(this->children != 0);
        if (this->isEndNode()) {
            this->deleteChildren();
        } else {
            for (int i = 0; i < getTDim(); i++) {
                assert(this->children[i] != 0);
                this->getChild(i).purgeGenerated();
            }
        }
    }
}

template<int D>
void MRNode<D>::getCenter(double *r) const {
    assert(r != 0);
    double sFac = pow(2.0, -getScale());
    for (int d = 0; d < D; d++) {
        double l = (double) getTranslation()[d];
        double o = this->tree->getOrigin()[d];
        r[d] = sFac*(l + 0.5) - o;
    }
}

template<int D>
void MRNode<D>::getLowerBounds(double *r) const {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRNode<D>::getUpperBounds(double *r) const {
    NOT_IMPLEMENTED_ABORT;
}

/** Given the child index (between 0 and (2^D - 1)) calculate the translation
  * index and store it in the *l argument. Assumes *l is allocated. */
template<int D>
void MRNode<D>::calcChildTranslation(int *transl, int cIdx) const {
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    const int *l = this->nodeIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        transl[d] = (2 * l[d]) + ((cIdx >> d) & 1);
    }
}

template<int D>
void MRNode<D>::calcChildNodeIndex(NodeIndex<D> &nIdx, int cIdx) const {
    int n = getScale() + 1;
    int rank = getRankId();
    int l[3];
    calcChildTranslation(l, cIdx);
    nIdx.setRankId(rank);
    nIdx.setScale(n);
    nIdx.setTranslation(l);
}

/** Routine to find the path along the tree.
  *
  * Given the translation indices at the final scale, computes the child m
  * to be followed at the current scale in oder to get to the requested
  * node at the final scale. The result is the index of the child needed.
  * The index is obtained by bit manipulation of of the translation indices. */
template<int D>
int MRNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    assert(isAncestor(nIdx));
    int cIdx = 0;
    int diffScale = nIdx.getScale() - getScale() - 1;
    assert(diffScale >= 0);
    for (int d = 0; d < D; d++) {
        int bit = (nIdx.getTranslation()[d] >> (diffScale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

/** Routine to find the path along the tree.
  *
  * Given a point in space, determines which child should be followed
  * to get to the corresponding terminal node. */
template<int D>
int MRNode<D>::getChildIndex(const double *r) const {
    assert(hasCoord(r));
    int cIdx = 0;
    double sFac = pow(2.0, -getScale());
    const double *o = this->tree->getOrigin();
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        if (r[d] > (sFac*(l[d] + 0.5) + o[d])) {
            cIdx = cIdx + (l[d] << d);
        }
    }
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    return cIdx;
}

/** Const version of node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
const MRNode<D> *MRNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
//    if (getScale() == idx.getScale()) { // we're done
//        return this;
//    }
//    if (this->isEndNode()) { // don't return GenNodes
//        return 0;
//    }
//    int cIdx = getChildIndex(idx);
//    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

/** Node retriever that NEVER generates.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or a NULL pointer if
  * the node does not exist, or if it is a GenNode. Recursion starts at at this
  * node and ASSUMES the requested node is in fact decending from this node. */
template<int D>
MRNode<D> *MRNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        return this;
    }
    if (this->isEndNode()) { // don't return GenNodes
        return 0;
    }
    int cIdx = getChildIndex(idx);
    return this->children[cIdx]->retrieveNodeNoGen(idx);
}

template<int D>
const MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const double *r, int depth) const {
    NOT_IMPLEMENTED_ABORT;
//    if (getDepth() == depth or this->isEndNode()) {
//        return this;
//    }
//    int cIdx = getChildIndex(r);
//    const MRNode<D> &child = getChild(cIdx);
//    return child.retrieveNodeOrEndNode(r, depth);
}

/** Node retriever that return requested ProjectedNode or EndNode.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine returns the appropriate ProjectedNode, or the EndNode on the
  * path to the requested node, and will never create or return GenNodes.
  * Recursion starts at at this node and ASSUMES the requested node is in fact
  * decending from this node. */
template<int D>
MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const double *r, int depth) {
    NOT_IMPLEMENTED_ABORT;
//    if (getDepth() == depth or this->isEndNode()) {
//        return this;
//    }
//    int cIdx = getChildIndex(r);
//    MRNode<D> &child = getChild(cIdx);
//    return child.retrieveNodeOrEndNode(r, depth);
}

template<int D>
const MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT;
//    if (nodeIndex.scale() == idx.getScale()) { // we're done
//        return this;
//    }
//    // We should in principle lock before read, but it makes things slower,
//    // and the EndNode status does not change (normally ;)
//    if (isEndNode()) {
//        return this;
//    }
//    int cIdx = getChildIndex(idx);
//    assert(children[cIdx] != 0);
//    return children[cIdx]->retrieveNodeOrEndNode(idx);
}

/** Node retriever that ALWAYS returns the requested node.
  *
  * Recursive routine to find and return the node with a given NodeIndex.
  * This routine always returns the appropriate node, and will generate nodes
  * that does not exist. Recursion starts at at this node and ASSUMES the
  * requested node is in fact decending from this node. */
template<int D>
MRNode<D> *MRNode<D>::retrieveNode(int n, const double *r) {
    NOT_IMPLEMENTED_ABORT;
//    if (this->nodeIndex.scale() == n) {
//        return this;
//    }
//    int l[D];
//    int idx = getChildIndex(r);
//    calcChildTranslation(idx, l);
//    // If we have reached an endNode, lock if necessary, and start generating
//    // NB! retrieveNode() for GenNodes behave a bit differently.
//    if (this->isEndNode()) {
//        SET_NODE_LOCK();
//        if (this->isLeafNode()) {
//            genChildren();
//        }
//        UNSET_NODE_LOCK();
//        return this->children[idx];
//    }
//    return this->children[idx]->retrieveNode(n, r);
}

/** Node retriever that ALWAYS returns the requested node, possibly without coefs.
  *
  * Recursive routine to find and return the node with a given NodeIndex. This
  * routine always returns the appropriate node, and will generate nodes that
  * does not exist. Recursion starts at at this node and ASSUMES the requested
  * node is in fact decending from this node. */
template<int D>
MRNode<D> *MRNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (getScale() == idx.getScale()) { // we're done
        return this;
    }
    if (isEndNode()) {
        SET_NODE_LOCK();
        if (isLeafNode()) {
            genChildren();
        }
        UNSET_NODE_LOCK();
    }
    int cIdx = getChildIndex(idx);
    return this->children[cIdx]->retrieveNode(idx);
}

/** Test if a given coordinate is within the boundaries of the node. */
template<int D>
bool MRNode<D>::hasCoord(const double *r) const {
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    const double *o = this->tree->getOrigin();
    println(1, "[" << r[0] << "," << r[1] << "," << r[2] << "]");
    println(1, "[" << o[0] << "," << o[1] << "," << o[2] << "]");
    println(1, "[" << l[0] << "," << l[1] << "," << l[2] << "]");
    println(1, *this);
    for (int d = 0; d < D; d++) {
        if (r[d] < (sFac*l[d] + o[d]) or r[d] > (sFac*(l[d] + 1)) + o[d]) {
            println(1, "false");
            return false;
        }
    }
    println(1, "true");
    return true;
}

/** Testing if nodes are compatible wrt NodeIndex and Tree (order, rootScale,
  * relPrec, etc). */
template<int D>
bool MRNode<D>::isCompatible(const MRNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
//    if (nodeIndex != node.nodeIndex) {
//        println(0, "nodeIndex mismatch" << std::endl);
//        return false;
//    }
//    if (not this->tree->checkCompatible(*node.tree)) {
//        println(0, "tree type mismatch" << std::endl);
//        return false;
//    }
//    return true;
}

/** Test if the node is decending from a given NodeIndex, that is, if they have
  * overlapping support. */
template<int D>
bool MRNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int relScale = idx.getScale() - getScale();
    if (relScale < 0) {
        return false;
    }
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx.getTranslation()[d] >> relScale;
        if (l[d] != reqTransl) {
            return false;
        }
    }
    return true;
}

template<int D>
bool MRNode<D>::isDecendant(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
void MRNode<D>::broadcastCoefs(int src, mpi::communicator *comm) {
    NOT_IMPLEMENTED_ABORT;
//#ifdef HAVE_MPI

//    if (comm != 0) {
//        comm = &node_group;
//    }
//    assert(this->isAllocated());
//    double *data = coefs->data();
//    mpi::broadcast(*comm, data, getNCoefs(), src);
//    this->setHasCoefs(true);
//    this->setFullRedundancy();
//#endif
//}

//template<int D>
//void MWNode<D>::sendCoefs(int dest, int tag,
//                             mpi::communicator *comm) {
//#ifdef HAVE_MPI

//    if (comm != 0) {
//        comm = &node_group;
//    }
//    const double *data = coefs->data();
//    comm->send(dest, tag, data, getNCoefs());
//    this->setRedundancy(dest);
//#endif
}

template<int D>
void MRNode<D>::assignDecendantTags(int rank) {
    for (int n = 0; n < getNChildren(); n++) {
        MRNode<D> &child = getChild(n);
        child.setRankId(rank);
        child.assignDecendantTags(rank);
    }
}

template class MRNode<1>;
template class MRNode<2>;
template class MRNode<3>;

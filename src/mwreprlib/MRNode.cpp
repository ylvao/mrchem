#include "MRNode.h"
#include "MathUtils.h"

using namespace std;

template<int D>
MRNode<D>::MRNode() {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>::MRNode(MRTree<D> *t, const NodeIndex<D> &idx) : nodeIndex(idx) {
    if (t == 0) {
        MSG_FATAL("Cannot initialize node without tree!");
    }

    this->tree = t;
    this->tree->incrementNodeCount(idx.getScale());

    this->parent = 0;
    this->status = 0;
    this->children = 0;

    setIsLeafNode();
    setIsEndNode();
    setIsRootNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
MRNode<D>::MRNode(MRNode<D> *p, const int *l) {
    this->parent = p;
    this->status = 0;
    this->children = 0;

    if (this->parent == 0) {
        NOT_IMPLEMENTED_ABORT;
    } else {
        int n = this->parent->getScale() + 1;
        this->nodeIndex.setScale(n);
        this->nodeIndex.setTranslation(l);
        this->tree = this->parent->tree;
        this->tree->incrementNodeCount(n);
    }
    setIsLeafNode();
#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
MRNode<D>::MRNode(const MRNode<D> &nd, MRNode<D> *p, bool copyCoefs) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D>::MRNode(const MRNode<D> &nd, MRTree<D> *t, bool copyCoefs) {
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
    NOT_IMPLEMENTED_ABORT;
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

template<int D>
void MRNode<D>::calcChildTranslation(int cIdx, int *transl) const {
    assert(cIdx >= 0);
    assert(cIdx < getTDim());
    const int *l = this->nodeIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        transl[d] = (2 * l[d]) + ((cIdx >> d) & 1);
    }
}

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

template<int D>
const MRNode<D> *MRNode<D>::retrieveNodeNoGen(const NodeIndex<D> &idx) const {
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
    if (getDepth() == depth or this->isEndNode()) {
        return this;
    }
    int cIdx = getChildIndex(r);
    const MRNode<D> &child = getChild(cIdx);
    return child.retrieveNodeOrEndNode(r, depth);
}

template<int D>
MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const double *r, int depth) {
    println(0, *this);
    if (getDepth() == depth or this->isEndNode()) {
        println(0, "returning");
        return this;
    }
    int cIdx = getChildIndex(r);
    println(0, "cIdx " << cIdx);
    MRNode<D> &child = getChild(cIdx);
    return child.retrieveNodeOrEndNode(r, depth);
}

template<int D>
const MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNode<D> *MRNode<D>::retrieveNodeOrEndNode(const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
bool MRNode<D>::hasCoord(const double *r) const {
    double sFac = pow(2.0, -getScale());
    const int *l = getTranslation();
    const double *o = this->tree->getOrigin();
    println(0, "[" << r[0] << "," << r[1] << "," << r[2] << "]");
    println(0, "[" << o[0] << "," << o[1] << "," << o[2] << "]");
    println(0, "[" << l[0] << "," << l[1] << "," << l[2] << "]");
    println(0, *this);
    for (int d = 0; d < D; d++) {
        if (r[d] < (sFac*l[d] + o[d]) or r[d] > (sFac*(l[d] + 1)) + o[d]) {
            println(0, "false");
            return false;
        }
    }
    println(0, "true");
    return true;
}

template<int D>
bool MRNode<D>::isCompatible(const MRNode<D> &node) {
    NOT_IMPLEMENTED_ABORT;
}

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
}

template<int D>
void MRNode<D>::sendCoefs(int who, int tag, mpi::communicator *comm) {
    NOT_IMPLEMENTED_ABORT;
}
template<int D>
void MRNode<D>::receiveCoefs(int who, int tag, mpi::communicator *comm) {

    NOT_IMPLEMENTED_ABORT;
}

template<int D>
mpi::request MRNode<D>::isendCoefs(int who, int tag, int comp) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
mpi::request MRNode<D>::ireceiveCoefs(int who, int tag, int comp) {
    NOT_IMPLEMENTED_ABORT;
}

template class MRNode<1>;
template class MRNode<2>;
template class MRNode<3>;

/*
template<int D>
GridNode<D> &GridNode<D>::getNode(const NodeIndex<D> &idx) {
    if (not isAncestor(idx)) {
 THROW_ERROR("Cannot get node, node out of bounds:\nthis->idx: " <<
                        this->getNodeIndex() << "\narg->idx: " << idx);
    }
    int scale = idx.getScale();
    if (scale > this->tree->getMaxScale()) {
 THROW_ERROR("Requested (" << scale << ") node beyond max scale ("
                        << this->tree->getMaxScale() << ")!")
    }
    return *(retrieveNode(idx));
}

template<int D>
GridNode<D> *GridNode<D>::retrieveNode(const NodeIndex<D> &idx) {
    if (this->nodeIndex.getScale() == idx.getScale()) { // we're done
        return this;
    }
    if (this->isLeafNode()) {
 THROW_ERROR("Requested node does not exist");
    }
    int childIdx = getChildIndex(idx);
    return children[childIdx]->retrieveNode(idx);
}

template<int D>
bool GridNode<D>::isAncestor(const NodeIndex<D> &idx) const {
    int inpScale = idx.getScale();
    int relScale = inpScale - this->nodeIndex.getScale();
    if (relScale < 0) {
        return false;
    }
    const int *l = this->nodeIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        int reqTransl = idx.getTranslation()[d] >> relScale;
 if (l[d] != reqTransl) {
     return false;
 }
    }
    return true;
}

template<int D>
int GridNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    int cIdx = 0;
    int delta_scale = nIdx.getScale() - this->nodeIndex.getScale() - 1;
    for (int d = 0; d < D; d++) {
        int bit = (nIdx.getTranslation()[d] >> (delta_scale)) & 1;
        cIdx = cIdx + (bit << d);
    }
    return cIdx;
}

template<int D>
NodeIndex<D> GridNode<D>::getChildIndex(int cIdx) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GridNode<D>::calcChildTranslation(int cIdx, int *transl) const {
    const int *l = this->nodeIndex.getTranslation();
    for (int d = 0; d < D; d++) {
        transl[d] = (2 * l[d]) + ((cIdx >> d) & 1);
    }
}

template<int D>
void GridNode<D>::getCenter(double *r) const {
    if (r == 0) {
 THROW_ERROR("Invalid argument");
    }
    const double *origin = this->tree->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
 r[d] = sFac*(1.0*l[d] + 0.5) - origin[d];
    }
}

template<int D>
void GridNode<D>::getLowerBounds(double *r) const {
    if (r == 0) {
 THROW_ERROR("Invalid argument");
    }
    const double *origin = this->tree->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
 r[d] = sFac*(1.0*l[d]) - origin[d];
    }
}

template<int D>
void GridNode<D>::getUpperBounds(double *r) const {
    if (r == 0) {
 THROW_ERROR("Invalid argument");
    }
    const double *origin = this->tree->getOrigin();
    int n = getScale();
    double sFac = pow(2.0, -n);
    const int *l = getTranslation();
    for (int d = 0; d < D; d++) {
 r[d] = sFac*(1.0*l[d] + 1.0) - origin[d];
    }
}

template<int D>
void GridNode<D>::allocKindergarten() {
    if (this->children == 0) {
 int nChildren = getTDim();
 this->children = new GridNode<D> *[nChildren];
 for (int n = 0; n < nChildren; n++) {
     this->children[n] = 0;
 }
    }
}


*/


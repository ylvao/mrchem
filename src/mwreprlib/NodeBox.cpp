/*
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "NodeBox.h"
#include "MRNode.h"

using namespace std;

template<int D>
NodeBox<D>::NodeBox(const NodeIndex<D> &idx, const int *nb)
        : BoundingBox<D>(idx, nb),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
NodeBox<D>::NodeBox(const BoundingBox<D> &box)
        : BoundingBox<D>(box),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
NodeBox<D>::NodeBox(const NodeBox<D> &box)
        : BoundingBox<D>(box),
          nOccupied(0),
          nodes(0) {
    allocNodePointers();
}

template<int D>
void NodeBox<D>::allocNodePointers() {
    assert(this->nodes == 0);
    int nNodes = this->size();
    this->nodes = new MRNode<D>*[nNodes];
    for (int n = 0; n < nNodes; n++) {
        this->nodes[n] = 0;
    }
    this->nOccupied = 0;
}

template<int D>
NodeBox<D>::~NodeBox() {
    deleteNodes();
}

template<int D>
void NodeBox<D>::deleteNodes() {
    if (this->nodes == 0) {
        return;
    }
    for (int n = 0; n < this->size(); n++) {
        removeNode(n);
    }
    delete [] this->nodes;
    this->nodes = 0;
}

template<int D>
void NodeBox<D>::setNode(int bIdx, MRNode<D> **node) {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    removeNode(bIdx);
    this->nodes[bIdx] = *node;
    this->nOccupied++;
    assert(this->nOccupied > 0);
    *node = 0;
}

/** Remove a node from the box **/
template<int D>
void NodeBox<D>::removeNode(int bIdx) {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    if (this->nodes[bIdx] != 0 ) {
        delete nodes[bIdx];
        this->nodes[bIdx] = 0;
        this->nOccupied--;
    }
    assert(this->nOccupied >= 0);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(const NodeIndex<D> &nIdx) {
    NOT_IMPLEMENTED_ABORT;
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(const double *r) {
    NOT_IMPLEMENTED_ABORT;
    int bIdx = this->getBoxIndex(r);
    return getNode(bIdx);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(int bIdx) {
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template<int D>
const MRNode<D>& NodeBox<D>::getNode(const NodeIndex<D> &nIdx) const {
    NOT_IMPLEMENTED_ABORT;
    int bIdx = this->getBoxIndex(nIdx);
    return getNode(bIdx);
}

template<int D>
const MRNode<D>& NodeBox<D>::getNode(const double *r) const {
    NOT_IMPLEMENTED_ABORT;
    int bIdx = this->getBoxIndex(r);
    return getNode(bIdx);
}

template<int D>
const MRNode<D>& NodeBox<D>::getNode(int bIdx) const {
    NOT_IMPLEMENTED_ABORT;
    assert(bIdx >= 0);
    assert(bIdx < this->nBoxes[D]);
    assert(this->nodes[bIdx] != 0);
    return *this->nodes[bIdx];
}

template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;

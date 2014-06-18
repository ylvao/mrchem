/*
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "NodeBox.h"
#include "MRNode.h"
#include "TelePrompter.h"
#include "MathUtils.h"
#include "constants.h"

using namespace std;

template<int D>
NodeBox<D>::NodeBox(int n, const int *nb, const double *o) 
    : BoundingBox<D>(n, nb, o), nodes(0) {
    NOT_IMPLEMENTED_ABORT
    allocNodePointers();
}

template<int D>
NodeBox<D>::NodeBox(const NodeBox<D> &box) 
	: BoundingBox<D>(box), nodes(0) {
    NOT_IMPLEMENTED_ABORT
    this->unitLength = box.unitLength;
    this->nBoxes[D] = box.nBoxes[D];
    for (int i = 0; i < D; i++) {
	assert(box.nBoxes[i] > 0);
	this->nBoxes[i] = box.nBoxes[i];
	this->origin[i] = box.origin[i];
	this->boxLength[i] = box.boxLength[i];
	this->lowerBounds[i] = box.lowerBounds[i];
	this->upperBounds[i] = box.upperBounds[i];
    }
    allocNodePointers();
}

template<int D>
void NodeBox<D>::allocNodePointers() {
    if(this->nodes != 0) {
	THROW_ERROR("Node pointers already allocated");
    }
    int n = this->getNBoxes(); 
    this->nodes = new MRNode<D>*[n];
    for (int i = 0; i < n; i++) {
	this->nodes[i] = 0;
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
    for (int i = 0; i < this->getNBoxes(); i++) {
	removeNode(i);
    }
    delete [] this->nodes;
    this->nodes = 0;
}

template<int D>
void NodeBox<D>::setNode(int idx, MRNode<D> **node) {
    if ((idx < 0) or (idx > this->nBoxes[D])) {
	MSG_FATAL("Node index out of range: " << idx);
    }
    removeNode(idx);
    this->nodes[idx] = *node;
    this->nOccupied++;
    *node = 0;
}

/** Remove a node from the box **/
template<int D>
void NodeBox<D>::removeNode(int idx) {
    assert(idx >= 0 and idx < this->nBoxes[D]);
    if (this->nodes[idx] != 0 ) {
	delete nodes[idx];
	this->nodes[idx] = 0;
	this->nOccupied--;
    }
    assert(this->nOccupied >= 0);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(const NodeIndex<D> &idx) {
    int i = getBoxIndex(idx);
    return getNode(i);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(const double *r) {
    int i = this->getBoxIndex(r);
    return getNode(i);
}

template<int D>
MRNode<D>& NodeBox<D>::getNode(int i) {
    assert(i >= 0 and i < this->nBoxes[D]);
    if (this->nodes[i] == 0) {
	println(1, *this);
	MSG_FATAL("Node not initialized, index: " << i);
    }
    return *this->nodes[i];
}

template class NodeBox<1>;
template class NodeBox<2>;
template class NodeBox<3>;

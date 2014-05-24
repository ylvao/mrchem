/*
 *  \date May 24, 2014
 *  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
 *          CTCC, University of Tromsø
 *
 */

#include "GridNodeBox.h"
#include "GridNode.h"
#include "TelePrompter.h"
#include "constants.h"

using namespace std;

template<int D>
GridNodeBox<D>::GridNodeBox(int scale, const int *nbox, const double *origo) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
GridNodeBox<D>::GridNodeBox(const GridNodeBox<D> &box) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GridNodeBox<D>::allocNodePointers() {
    if(this->nodes != 0) {
	THROW_ERROR("Node pointers already allocated");
    }
    int n = getNBoxes(); 
    this->nodes = new GridNode<D>*[n];
    for (int i = 0; i < n; i++) {
	this->nodes[i] = 0;
    }
    this->nOccupied = 0;
}

template<int D>
GridNodeBox<D>::~GridNodeBox() {
    deleteNodes();
}

template<int D>
void GridNodeBox<D>::deleteNodes() {
    NOT_IMPLEMENTED_ABORT
/*
    if (this->nodes == 0) {
	return;
    }
    for (int i = 0; i < this->nBoxes; i++) {
	if (this->nodes[i] != 0) {
	    delete this->nodes[i];
	}
	this->nodes[i] = 0;
    }
    delete [] this->nodes;
    this->nodes = 0;
*/
}

template<int D>
void GridNodeBox<D>::setNode(int idx, GridNode<D> *node) {
    NOT_IMPLEMENTED_ABORT
/*
	if ((idx < 0) or (idx > this->nBoxes)) {
		MSG_FATAL("Node index out of range: " << idx);
	}
	assert(nodes[idx] == 0);
	nodes[idx] = node;
	nOccupied++;
*/
}

/** Remove a node from the box **/
template<int D>
void GridNodeBox<D>::removeNode(int idx) {
    NOT_IMPLEMENTED_ABORT
/*
	assert(idx >= 0 and idx < this->nBoxes);
	assert(nodes[idx] != 0);
	if (tree != 0) {
		delete nodes[idx];
	}
	nodes[idx] = 0;
	nOccupied--;
*/
}

template<int D>
NodeIndex<D> GridNodeBox<D>::getNodeIndex(const double *r) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
NodeIndex<D> GridNodeBox<D>::getNodeIndex(int i) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int GridNodeBox<D>::getBoxIndex(const double *r) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int GridNodeBox<D>::getBoxIndex(const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT
}
    

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(const NodeIndex<D> &idx) const {
    NOT_IMPLEMENTED_ABORT
/*
	int rootIdx = this->getBoxIndex(idx);
	MWNode<D> &rootnode = getNode(rootIdx);
	return rootnode.getNode(idx, genEmpty);
*/
}

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(const double *r) const {
    NOT_IMPLEMENTED_ABORT
/*
    int i = this->getBoxIndex(r);
    return getNode(i);
*/
}

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(int i) const {
    NOT_IMPLEMENTED_ABORT
/*
    assert(i >= 0 and i < this->nBoxes);
    if (this->nodes[i] == 0) {
	println(1, *this);
	MSG_FATAL("Node not initialized, index: " << i);
    }
    return *this->nodes[i];
*/
}

template<int D>
int GridNodeBox<D>::getNBoxes(int d) const {
    if (d < 0) {
	return this->nBoxes[D];
    } else if (d < D) {
	return this->nBoxes[d];
    } else {
	THROW_ERROR("Invalid dimension argument");
    }
}

template class GridNodeBox<1>;
template class GridNodeBox<2>;
template class GridNodeBox<3>;

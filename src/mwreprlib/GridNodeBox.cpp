/*
 *  \date Aug 18, 2009
 *  \author Jonas Juselius <jonas.juselius@uit.no> \n
 *          CTCC, University of Tromsø
 *
 * \breif NodeBox provides a means to have multiple root nodes,
 * and thus select other shapes of the world than square, cubic or
 * hypercubic.
 *
 * NOTE! IMPORTANT! The box is created with N empty slots for
 * root nodes. If you change the box AFTER adding root nodes, all nodes
 * get deleted and the box has to be initialized with new,
 * fresh-from-the-factory root nodes. Batteries not included.
 */

#include <cmath>
#include <string>

#include "GridNodeBox.h"
#include "GridNode.h"
#include "TelePrompter.h"
#include "constants.h"

using namespace std;

/** The default constructor just initializes the box with copies the
 * default values as provided by the global (static) worldbox */
template<int D>
GridNodeBox<D>::GridNodeBox(const BoundingBox<D> *box) : BoundingBox<D>() {
    this->nodes = 0;
    this->nOccupied = 0;
    if (box == 0) {
	NOT_IMPLEMENTED_ABORT
    } else {
	copyBox(*box);
    }
    allocNodePointers(this->nBoxes);
}

template<int D>
GridNodeBox<D>::~GridNodeBox() {
    deleteNodes();
}

template<int D>
void GridNodeBox<D>::deleteNodes() {
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
    int i = this->getBoxIndex(r);
    return getNode(i);
}

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(int i) const {
    assert(i >= 0 and i < this->nBoxes);
    if (this->nodes[i] == 0) {
	println(1, *this);
	MSG_FATAL("Node not initialized, index: " << i);
    }
    return *this->nodes[i];
}

template<int D>
void GridNodeBox<D>::allocNodePointers(int n) {
    if(this->nodes != 0) {
	THROW_ERROR("Node pointers already allocated");
    }
    this->nodes = new GridNode<D>*[n];
    for (int i = 0; i < this->nBoxes; i++) {
	this->nodes[i] = 0;
    }
    this->nOccupied = 0;
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

template class GridNodeBox<1>;
template class GridNodeBox<2>;
template class GridNodeBox<3>;

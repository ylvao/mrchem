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
GridNodeBox<D>::GridNodeBox(int scale, const int *nbox, const double *origo) 
    : cornerIndex(scale), nodes(0) {
    initBox(nbox, origo);
    allocNodePointers();
}

template<int D>
GridNodeBox<D>::GridNodeBox(const GridNodeBox<D> &box) 
	: cornerIndex(box.cornerIndex), nodes(0) {
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
void GridNodeBox<D>::initBox(const int *nbox, const double *origo) {
    this->unitLength = 1.0/pow(2.0, this->cornerIndex.getScale());
    this->nBoxes[D] = 1;
    for (int i = 0; i < D; i++) {
	if (nbox == 0) {
	    this->nBoxes[i] = 1;
	} else {
	    assert(nbox[i] > 0);
	    this->nBoxes[i] = nbox[i];
        }
        if (origo == 0) {
	    this->origin[i] = 0.0;
        } else {
	    this->origin[i] = origo[i];
	}
	this->boxLength[i] = this->unitLength * this->nBoxes[i];
	this->nBoxes[D] *= this->nBoxes[i];
	this->lowerBounds[i] = 0.0 - this->origin[i];
	this->upperBounds[i] = this->boxLength[i] - this->origin[i];
    }
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
    if (this->nodes == 0) {
	return;
    }
    for (int i = 0; i < getNBoxes(); i++) {
	removeNode(i);
    }
    delete [] this->nodes;
    this->nodes = 0;
}

template<int D>
void GridNodeBox<D>::setNode(int idx, GridNode<D> **node) {
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
void GridNodeBox<D>::removeNode(int idx) {
    assert(idx >= 0 and idx < this->nBoxes[D]);
    if (this->nodes[idx] != 0 ) {
	delete nodes[idx];
	this->nodes[idx] = 0;
	this->nOccupied--;
    }
    assert(this->nOccupied >= 0);
}

template<int D>
NodeIndex<D> GridNodeBox<D>::getNodeIndex(const double *r) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
NodeIndex<D> GridNodeBox<D>::getNodeIndex(int bIdx) {
    assert(bIdx >= 0 and bIdx <= nBoxes[D]);
    NodeIndex<D> nIdx(getRootScale());
    if (D == 1) {
        int l = bIdx - this->cornerIndex.getTranslation()[0];
        nIdx.setTranslation(&l);
	return nIdx;
    }
    int transl[D];
    for (int i = D - 1; i >= 0; i--) {
        int ncells = 1;
        for (int j = 0; j < i; j++) {
            ncells *= this->nBoxes[j];
        }
        double div = bIdx / ncells;
        double iint;
        modf(div, &iint);
        transl[i] = (int) iint;
        bIdx -= ncells * transl[i];
    }
    for (int i = 0; i < D; i++) {
        transl[i] += this->cornerIndex[i];
    }
    nIdx.setTranslation(transl);
    return nIdx;
}

template<int D>
int GridNodeBox<D>::getBoxIndex(const double *r) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
int GridNodeBox<D>::getBoxIndex(const NodeIndex<D> &nIdx) {
    int bIdx = 0;
    int scale = nIdx.getScale();
    const int *l = nIdx.getTranslation();
    const int *corner_l = this->cornerIndex.getTranslation();
    assert(scale >= this->cornerIndex.getScale());
    int relScale = scale - this->cornerIndex.getScale();
    if (D == 1) {
        int reqTransl = (l[0] >> relScale) - corner_l[0];
        return reqTransl;
    }
    for (int i = D - 1; i >= 0; i--) {
        int ncells = 1;
        for (int j = 0; j < i; j++) {
            ncells *= this->nBoxes[j];
        }
        int reqTransl = (l[i] >> relScale) - corner_l[i];
        bIdx += ncells * reqTransl;
    }
    return bIdx;
}
    

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(const NodeIndex<D> &idx) {
    int rootIdx = getBoxIndex(idx);
    GridNode<D> &rootnode = getNode(rootIdx);
    return rootnode.getNode(idx);
}

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(const double *r) {
    NOT_IMPLEMENTED_ABORT
/*
    int i = this->getBoxIndex(r);
    return getNode(i);
*/
}

template<int D>
GridNode<D> &GridNodeBox<D>::getNode(int i) {
    assert(i >= 0 and i < this->nBoxes[D]);
    if (this->nodes[i] == 0) {
	println(1, *this);
	MSG_FATAL("Node not initialized, index: " << i);
    }
    return *this->nodes[i];
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

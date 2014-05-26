/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include "GridNode.h"
#include "NodeIndex.h"

template<int D>
GridNode<D>::GridNode(MRGrid<D> *_grid, int n, const int *l) {
    if (_grid == 0) {
        MSG_FATAL("Cannot initialize node without tree!");
    }
    this->nodeIndex.setScale(n);
    this->nodeIndex.setTranslation(l);

    this->grid = _grid;
    this->grid->incrementNodeCount(n);

    this->parent = 0;
    this->status = 0;
    this->children = 0;

    this->roots = 0;
    this->weights = 0;

    setIsLeafNode();
    setIsEndNode();
    setIsRootNode();

#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
GridNode<D>::GridNode(GridNode<D> *_parent, const int *l) {
    this->parent = _parent;
    this->status = 0;
    this->children = 0;
    this->roots = 0;
    this->weights = 0;

    if (this->parent == 0) {
	NOT_IMPLEMENTED_ABORT
    } else {
	int n = this->parent->getScale() + 1;
        this->nodeIndex.setScale(n);
	this->nodeIndex.setTranslation(l);
        this->grid = this->parent->grid;
	this->grid->incrementNodeCount(n);
    }
    setIsLeafNode();
    setIsEndNode();
#ifdef OPENMP
    omp_init_lock(&node_lock);
#endif
}

template<int D>
GridNode<D>::~GridNode() {
    SET_NODE_LOCK();
    if (this->roots != 0) {
        delete this->roots;
    }
    if (this->weights != 0) {
        delete this->weights;
    }
    if (this->isBranchNode()) {
        assert(this->children != 0);
        deleteChildren();
    }
    this->grid->decrementNodeCount(getScale());
    UNSET_NODE_LOCK();
#ifdef OPENMP
    omp_destroy_lock(&node_lock);
#endif
}

template<int D>
void GridNode<D>::deleteChildren() {
    assert(this->children != 0);
    for (int i = 0; i < this->tDim; i++) {
        if (this->children[i] != 0) {
            delete this->children[i];
            this->children[i] = 0;
        }
    }
    delete [] this->children;
    this->children = 0;
    this->setIsLeafNode();
}


template<int D>
void GridNode<D>::createChildren() {
    if (this->children == 0) {
        this->allocKindergarten();
    }
    for (int i = 0; i < this->tDim; i++) {
        createChild(i);
    }
    this->setIsBranchNode();
    this->clearIsEndNode();
}

template<int D>
void GridNode<D>::createChild(int i) {
    assert(this->children[i] == 0);
    int l[3];
    calcChildTranslation(i, l);
    GridNode<D> *child = new GridNode<D> (this, l);
    this->children[i] = child;
    child->setIsEndNode();
}

template<int D>
void GridNode<D>::allocKindergarten() {
    if (this->children == 0) {
	this->children = new GridNode<D> *[this->tDim];
	for (int i = 0; i < this->tDim; i++) {
	    this->children[i] = 0;
	}
    }
}

template<int D>
GridNode<D> &GridNode<D>::getNode(const NodeIndex<D> &idx) {
    if (not isAncestor(idx)) {
	THROW_ERROR("Cannot get node, node out of bounds:\nthis->idx: " <<
                        this->getNodeIndex() << "\narg->idx: " << idx);
    }
    int scale = idx.scale();
    if (scale > this->grid->getMaxScale()) {
	THROW_ERROR("Requested (" << scale << ") node beyond max scale ("
                        << this->grid->getMaxScale() << ")!")
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
    int inpScale = idx.scale();
    int relScale = inpScale - this->nodeIndex.scale();
    if (relScale < 0) {
        return false;
    }
    for (int i = 0; i < D; i++) {
        NodeIndex<D> &l = const_cast<NodeIndex<D> &> (idx);
        int reqTransl = l[i] >> relScale;
	if (this->nodeIndex[i] != reqTransl) {
	    return false;
	}
    }
    return true;
}

template<int D>
int GridNode<D>::getChildIndex(const NodeIndex<D> &nIdx) const {
    int cIdx = 0;
    int delta_scale = nIdx.getScale() - this->nodeIndex.getScale() - 1;
    for (int i = 0; i < D; i++) {
        int bit = (nIdx.getTranslation()[i] >> (delta_scale)) & 1;
        cIdx = cIdx + (bit << i);
    }
    return cIdx;
}

template<int D>
NodeIndex<D> GridNode<D>::getChildIndex(int cIdx) const {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GridNode<D>::calcChildTranslation(int cIdx, int *l) const {
    for (int i = 0; i < D; i++) {
        l[i] = (2 * this->nodeIndex[i]) + ((cIdx >> i) & 1);
    }
}

template class GridNode<1>;
template class GridNode<2>;
template class GridNode<3>;

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
GridNode<D>::GridNode(MRGrid<D> *_grid, const NodeIndex<D> &idx) {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
GridNode<D>::~GridNode() {
    NOT_IMPLEMENTED_ABORT
}

template<int D>
void GridNode<D>::deleteChildren() {
    NOT_IMPLEMENTED_ABORT
}


template<int D>
void GridNode<D>::createChildren() {
    NOT_IMPLEMENTED_ABORT
}

template class GridNode<1>;
template class GridNode<2>;
template class GridNode<3>;

/**
*
*
*  \date May 23, 2014
*  \author Stig Rune Jensen <stig.r.jensen@uit.no> \n
*   CTCC, University of Tromsø
*
*/

#include "GridGenerator.h"
#include "TelePrompter.h"

template<int D>
GridGenerator<D>::~GridGenerator() {
    if (this->grid != 0) {
        THROW_ERROR("Grid pointer not released");
    }
}

template<int D>
void GridGenerator<D>::generateGrid(MRGrid<D> &outGrid) { 
    this->grid = &outGrid;
    NOT_IMPLEMENTED_ABORT;
    this->grid = 0;
}

template<int D>
void GridGenerator<D>::splitNodeTable(GridNodeVector &nodeTable) { 
    NOT_IMPLEMENTED_ABORT
}

template<int D>
bool GridGenerator<D>::splitCheck(GridNode<D> *node) {
    NOT_IMPLEMENTED_ABORT
}

template class GridGenerator<1>;
template class GridGenerator<2>;
template class GridGenerator<3>;

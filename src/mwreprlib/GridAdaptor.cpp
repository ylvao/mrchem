#include "GridAdaptor.h"
#include "MRGrid.h"
#include "FunctionTree.h"

template<int D>
GridAdaptor<D>::GridAdaptor(double pr, bool abs) {
    this->prec = pr;
    this->absPrec = abs;
}

template<int D>
void GridAdaptor<D>::adaptGrid(MRGrid<D> &grid, FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
MRNodeVector& GridAdaptor<D>::splitNodeTable(MRNodeVector &nodeTable) {
    NOT_IMPLEMENTED_ABORT;
}

template class GridAdaptor<1>;
template class GridAdaptor<2>;
template class GridAdaptor<3>;

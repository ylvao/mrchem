#include "GridAdaptor.h"
#include "MRGrid.h"
#include "MRNode.h"
#include "FunctionTree.h"

template<int D>
GridAdaptor<D>::GridAdaptor(double pr, bool abs) : prec(pr), absPrec(abs) {
}

template<int D>
void GridAdaptor<D>::adaptGrid(MRGrid<D> &outGrid, FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template<int D>
bool GridAdaptor<D>::splitCheck(const MRNode<D> *node) {
    NOT_IMPLEMENTED_ABORT;
}

template class GridAdaptor<1>;
template class GridAdaptor<2>;
template class GridAdaptor<3>;

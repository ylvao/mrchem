#include "FunctionProjector.h"
#include "RepresentableFunction.h"
#include "FunctionTree.h"
#include "GridAdaptor.h"
#include "MRGrid.h"

template<int D>
FunctionProjector<D>::FunctionProjector(MRGrid<D> *startGrid,
        GridAdaptor<D> *adap) : MWProjector<D>(startGrid, adap) {
}

template<int D>
void FunctionProjector<D>::operator()(FunctionTree<D> &tree,
                                   RepresentableFunction<D> &func) {
    NOT_IMPLEMENTED_ABORT;
}

template class FunctionProjector<1>;
template class FunctionProjector<2>;
template class FunctionProjector<3>;

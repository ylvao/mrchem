#include "MWProjector.h"
#include "GridAdaptor.h"
#include "MRGrid.h"
#include "FunctionTree.h"

template<int D>
MWProjector<D>::MWProjector() {
    this->adaptor = 0;
}

template<int D>
MWProjector<D>::MWProjector(GridAdaptor<D> &a) {
    this->adaptor = &a;
}

template<int D>
MWProjector<D>::~MWProjector() {
    this->adaptor = 0;
}

template<int D>
void MWProjector<D>::buildTree(FunctionTree<D> &tree) {
    NOT_IMPLEMENTED_ABORT;
}

template class MWProjector<1>;
template class MWProjector<2>;
template class MWProjector<3>;

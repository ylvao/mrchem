#include "FunctionProjector.h"
#include "RepresentableFunction.h"
#include "FunctionTree.h"
#include "GridAdaptor.h"
#include "MRGrid.h"

template<int D>
FunctionProjector<D>::FunctionProjector() {
    this->function = 0;
}

template<int D>
FunctionProjector<D>::FunctionProjector(GridAdaptor<D> &a) : MWProjector<D>(a) {
    this->function = 0;
}

template<int D>
FunctionProjector<D>::~FunctionProjector() {
    if (this->function != 0) {
        MSG_ERROR("Projector not properly cleared");
    }
}

template<int D>
void FunctionProjector<D>::operator()(FunctionTree<D> &tree,
                                      RepresentableFunction<D> &func) {
    this->function = &func;
    this->buildTree(tree);
    tree.mwTransformUp();
    this->function = 0;
}

template class FunctionProjector<1>;
template class FunctionProjector<2>;
template class FunctionProjector<3>;

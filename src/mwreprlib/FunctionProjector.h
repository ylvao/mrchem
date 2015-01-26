#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "MWProjector.h"

template<int D> class FunctionTree;
template<int D> class RepresentableFunction;

template<int D>
class FunctionProjector : public MWProjector<D> {
public:
    FunctionProjector(MRGrid<D> *startGrid = 0, GridAdaptor<D> *adap = 0);

    void operator()(FunctionTree<D> &tree, RepresentableFunction<D> &func);
protected:
};

#endif // FUNCTIONPROJECTOR_H

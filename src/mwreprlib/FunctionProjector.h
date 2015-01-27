#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "MWProjector.h"

template<int D> class FunctionTree;
template<int D> class RepresentableFunction;

template<int D>
class FunctionProjector : public MWProjector<D> {
public:
    FunctionProjector();
    FunctionProjector(GridAdaptor<D> &a);
    virtual ~FunctionProjector();

    void operator()(FunctionTree<D> &tree, RepresentableFunction<D> &func);

protected:
    RepresentableFunction<D> *function;

};

#endif // FUNCTIONPROJECTOR_H

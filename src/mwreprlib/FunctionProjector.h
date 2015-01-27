#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "MWProjector.h"
#include "mwrepr_declarations.h"

template<int D>
class FunctionProjector : public MWProjector<D> {
public:
    FunctionProjector();
    FunctionProjector(GridAdaptor<D> &a);
    virtual ~FunctionProjector();

    void operator()(FunctionTree<D> &t, RepresentableFunction<D> &f);

protected:
    RepresentableFunction<D> *func;

    void calcWaveletCoefs(MWNode<D> &node);
};

#endif // FUNCTIONPROJECTOR_H

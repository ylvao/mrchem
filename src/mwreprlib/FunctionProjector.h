#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "MWProjector.h"
#include "mwrepr_declarations.h"

template<int D>
class FunctionProjector : public MWProjector<D> {
public:
    FunctionProjector();
    FunctionProjector(MWAdaptor<D> &a);
    virtual ~FunctionProjector();

    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp);

protected:
    RepresentableFunction<D> *inpFunc;

    void calcWaveletCoefs(MWNode<D> &node);
};

#endif // FUNCTIONPROJECTOR_H

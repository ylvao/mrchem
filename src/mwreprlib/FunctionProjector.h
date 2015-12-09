#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "TreeBuilder.h"
#include "AnalyticProjector.h"
#include "FunctionTree.h"

template<int D>
class FunctionProjector : public TreeBuilder<D> {
public:
    FunctionProjector(TreeAdaptor<D> &a, int iter = -1)
        : TreeBuilder<D>(a, iter) { }
    virtual ~FunctionProjector() { }

    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
        this->projector = new AnalyticProjector<D>(inp);
        this->build(out);
        out.mwTransformUp();
        this->clearProjector();
    }
};

#endif // FUNCTIONPROJECTOR_H

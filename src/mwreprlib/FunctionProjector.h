#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "TreeBuilder.h"
#include "AnalyticCalculator.h"
#include "FunctionTree.h"

template<int D>
class FunctionProjector : public TreeBuilder<D> {
public:
    FunctionProjector(TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(iter) {
        NOT_IMPLEMENTED_ABORT;
//        this->adaptor = new TreeAdaptor<D>(a);
    }
    virtual ~FunctionProjector() {
        NOT_IMPLEMENTED_ABORT;
        this->clearAdaptor();
    }

    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
        NOT_IMPLEMENTED_ABORT;
        this->calculator = new AnalyticCalculator<D>(inp);
        this->build(out);
        out.mwTransform(BottomUp);
        this->clearCalculator();
    }
};

#endif // FUNCTIONPROJECTOR_H

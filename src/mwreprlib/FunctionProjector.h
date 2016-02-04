#ifndef FUNCTIONPROJECTOR_H
#define FUNCTIONPROJECTOR_H

#include "TreeBuilder.h"
#include "AnalyticCalculator.h"
#include "FunctionTree.h"

template<int D>
class FunctionProjector : public TreeBuilder<D> {
public:
    FunctionProjector()
            : TreeBuilder<D>(-1) {
        this->adaptor = new TreeAdaptor<D>();
    }
    FunctionProjector(const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(iter) {
       this->adaptor = a.copy();
    }
    virtual ~FunctionProjector() {
        this->clearAdaptor();
    }

    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
        this->calculator = new AnalyticCalculator<D>(inp);
        this->build(out);
        out.mwTransform(BottomUp);
        this->clearCalculator();
    }
};

#endif // FUNCTIONPROJECTOR_H

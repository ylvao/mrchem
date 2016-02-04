#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator(int iter = -1)
            : TreeBuilder<D>(iter) {
        this->calculator = new TreeCalculator<D>();
    }
    virtual ~GridGenerator() {
        this->clearCalculator();
    }

//    void operator()(FunctionTree<D> &out, RepresentableFunction<D> &inp) {
//        this->adaptor = new AnalyticAdaptor<D>(inp);
//        this->clearAdaptor();
//        NOT_IMPLEMENTED_ABORT;
//    }
    void operator()(FunctionTree<D> &out, MRTree<D> &inp) {
//        this->adaptor = new CopyAdaptor<D>(inp);
        NOT_IMPLEMENTED_ABORT;
        this->clearAdaptor();
    }
};

#endif // GRIDGENERATOR_H

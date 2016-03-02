#ifndef GRIDGENERATOR_H
#define GRIDGENERATOR_H

#include "TreeBuilder.h"

template<int D>
class GridGenerator : public TreeBuilder<D> {
public:
    GridGenerator(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->calculator = new TreeCalculator<D>();
    }
    virtual ~GridGenerator() {
        this->clearCalculator();
    }

    FunctionTree<D> *operator()() {
        return new FunctionTree<D>(this->MRA);
    }

    FunctionTree<D> *operator()(const RepresentableFunction<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    FunctionTree<D> *operator()(const FunctionTree<D> &inp) {
        FunctionTree<D> *out = new FunctionTree<D>(this->MRA);
        (*this)(*out, inp);
        return out;
    }

    void operator()(FunctionTree<D> &out, const RepresentableFunction<D> &inp) {
        NOT_IMPLEMENTED_ABORT;
//        this->adaptor = new AnalyticAdaptor<D>(inp);
//        this->clearAdaptor();
    }

    void operator()(FunctionTree<D> &out, const FunctionTree<D> &inp) {
        NOT_IMPLEMENTED_ABORT;
//        this->adaptor = new CopyAdaptor<D>(inp);
//        this->clearAdaptor();
    }
};

#endif // GRIDGENERATOR_H

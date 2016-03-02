#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "TreeBuilder.h"

template<int D>
class GridCleaner : public TreeBuilder<D> {
public:
    GridCleaner(const MultiResolutionAnalysis<D> &mra, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->calculator = new TreeCalculator<D>();
        this->adaptor = new TreeAdaptor<D>();
    }
    GridCleaner(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a, int iter = -1)
            : TreeBuilder<D>(mra, iter) {
        this->calculator = new TreeCalculator<D>();
        this->adaptor = a.copy();
    }
    virtual ~GridCleaner() {
        this->clearCalculator();
    }

    int operator()(FunctionTree<D> &out) {
        NOT_IMPLEMENTED_ABORT;
    }
};

#endif // GRIDCLEANER_H

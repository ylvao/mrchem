#ifndef GRIDCLEANER_H
#define GRIDCLEANER_H

#include "TreeBuilder.h"

template<int D>
class GridCleaner : public TreeBuilder<D> {
public:
    GridCleaner(const MultiResolutionAnalysis<D> &mra)
            : TreeBuilder<D>(mra, -1) {
        this->calculator = new TreeCalculator<D>();
    }
    GridCleaner(const MultiResolutionAnalysis<D> &mra,
                const TreeAdaptor<D> &a)
            : TreeBuilder<D>(mra, -1) {
        this->calculator = new TreeCalculator<D>();
        this->adaptor = a.copy();
    }
    virtual ~GridCleaner() {
        this->clearCalculator();
        this->clearAdaptor();
    }

    int operator()(FunctionTree<D> &out) {
        return clean(out);
    }

protected:
    int clean(MWTree<D> &tree) {
        if (this->calculator == 0) MSG_ERROR("Calculator not initialized");
        println(10, " == Clearing tree");

        int nSplit = 0;
        if (this->adaptor != 0) {
            MWNodeVector *endVec = tree.copyEndNodeTable();
            endVec = this->clearForeignNodes(endVec);
            MWNodeVector *splitVec = this->adaptor->splitNodeVector(*endVec);
            NodeIndexSet *splitSet = this->getNodeIndexSet(*splitVec);
            nSplit = splitVec->size();

            broadcast_index_list<D>(*splitSet);
            tree.splitNodes(*splitSet);//allocate new nodes

            delete splitVec;
            delete splitSet;
        }

        printout(10, "  -- #  0: Split        " << std::setw(6) << nSplit << " nodes\n");
        printout(10, "  -- #  1: Cleared      ");

        MWNodeVector nodeVec;
        tree.makeNodeTable(nodeVec);
        this->calculator->calcNodeVector(nodeVec);//clear all coefficients

        tree.resetEndNodeTable();
        tree.clearSquareNorm();
        return nSplit;
    }
};

#endif // GRIDCLEANER_H

#ifndef TREEBUILDER_H
#define TREEBUILDER_H

#include "mwrepr_declarations.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder(int iter);
    virtual ~TreeBuilder();

    void setMaxIter(int iter) { this->maxIter = iter; }

protected:
    int maxIter;
    TreeAdaptor<D> *adaptor;
    TreeCalculator<D> *calculator;

    void clearCalculator();
    void clearAdaptor();

    void build(MWTree<D> &tree);

    MWNodeVector* clearForeignNodes(MWNodeVector *oldVec) const;
    NodeIndexSet* getNodeIndexSet(const MWNodeVector &nodeVec) const;

    bool maxIterReached(int iter) const {
        if (this->maxIter < 0) return false;
        if (this->maxIter > iter) return false;
        return true;
    }
};

#endif // TREEBUILDER_H

#ifndef TREEBUILDER_H
#define TREEBUILDER_H

#include "mwrepr_declarations.h"
#include "TreeAdaptor.h"
#include "TreeProjector.h"

template<int D>
class TreeBuilder {
public:
    TreeBuilder(TreeAdaptor<D> &a, int iter);
    virtual ~TreeBuilder();

    void setMaxIter(int iter) { this->maxIter = iter; }

protected:
    int maxIter;
    TreeAdaptor<D> *adaptor;
    TreeProjector<D> *projector;

    void clearProjector();
    void build(MWTree<D> &tree);

    MRNodeVector* clearForeignNodes(MRNodeVector *oldVec) const;
    NodeIndexSet* getNodeIndexSet(const MRNodeVector &nodeVec) const;

    bool maxIterReached(int iter) const {
        if (this->maxIter < 0) return false;
        if (this->maxIter > iter) return false;
        return true;
    }
};

#endif // TREEBUILDER_H

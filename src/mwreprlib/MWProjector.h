#ifndef MWPROJECTOR_H
#define MWPROJECTOR_H

#include "mwrepr_declarations.h"
#include "MWAdaptor.h"

template<int D>
class MWProjector {
public:
    MWProjector(const MWAdaptor<D> &a, int iter) : adaptor(a), maxIter(iter) { }
    virtual ~MWProjector() { }

    void setMaxIter(int iter) { this->maxIter = iter; }
    MWAdaptor<D> &getAdaptor() { return this->adaptor; }

protected:
    int maxIter;
    MWAdaptor<D> adaptor;

    void buildTree(MWTree<D> &outTree);
    void calcNodeVector(MRNodeVector &nodeVec);
    virtual void calcNode(MWNode<D> &node) = 0;

    MRNodeVector* clearForeignNodes(MRNodeVector *oldVec) const;
    NodeIndexSet* getNodeIndexSet(const MRNodeVector &nodeVec) const;

    bool maxIterReached(int iter) const {
        if (this->maxIter < 0) return false;
        if (this->maxIter > iter) return false;
        return true;
    }
};

#endif // MWPROJECTOR_H
